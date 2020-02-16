/*
Copyright (C) 2019, NVIDIA CORPORATION. All rights reserved.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/


// Testing of panini... replace with fov cvar
#define fovx 180




// Forward transforms are used to calculate a known ray direction to a screen position (only used for temporal denoising etc.)
// Reverse transforms are used to get a ray direction from a screen position (used for the main path tracing)

void view_to_latlon(vec3 view, out float lat, out float lon)
{
   lon = atan(view.x, view.z);
   lat = atan(-view.y, sqrt(view.x*view.x+view.z*view.z));
}

void latlon_to_view(float lat, float lon, out vec3 view)
{
   float clat = cos(lat);
   view.x = sin(lon)*clat;
   view.y = -sin(lat);
   view.z = cos(lon)*clat;
}
							   

// ---------- Rectilinear transforms ----------
												  

bool rectilinear_forward(vec3 view_pos, out vec2 screen_pos, out float distance, bool previous)
{
	vec4 clip_pos;
	if(previous)
		clip_pos = global_ubo.P_prev * vec4(view_pos, 1);
	else
		clip_pos = global_ubo.P * vec4(view_pos, 1);

	vec3 normalized = clip_pos.xyz / clip_pos.w;
	screen_pos.xy = normalized.xy * 0.5 + vec2(0.5);
	distance = length(view_pos);

	return screen_pos.y > 0 && screen_pos.y < 1 && screen_pos.x > 0 && screen_pos.x < 1 && view_pos.z > 0;
}

vec3 rectlinear_reverse(vec2 screen_pos, float distance, bool previous)
{
	vec4 clip_pos = vec4(screen_pos.xy * 2.0 - vec2(1.0), 1, 1);
	vec3 view_dir;
	if(previous)
		view_dir = normalize((global_ubo.invP_prev * clip_pos).xyz);
	else
		view_dir = normalize((global_ubo.invP * clip_pos).xyz);
				
	return view_dir * distance;
}

// ---------- Cylindrical transforms ----------
bool cylindrical_forward(vec3 view_pos, out vec2 screen_pos, out float distance, bool previous, float cylindrical_hfov)
{
	float y = view_pos.y / length(view_pos.xz);
	if(previous)
		y *= global_ubo.P_prev[1][1];
	else
		y *= global_ubo.P[1][1];
	screen_pos.y = y * 0.5 + 0.5;

	float angle = atan(view_pos.x, view_pos.z);
	screen_pos.x = (angle / cylindrical_hfov) + 0.5;

	distance = length(view_pos);

	return screen_pos.y > 0 && screen_pos.y < 1 && screen_pos.x > 0 && screen_pos.x < 1;
}

vec3 cylindrical_reverse(vec2 screen_pos, float distance, bool previous, float cylindrical_hfov)
{
	vec4 clip_pos = vec4(0, screen_pos.y * 2.0 - 1.0, 1, 1);
	vec3 view_dir;
	if(previous)
		view_dir.y = (global_ubo.invP_prev * clip_pos).y;
	else
		view_dir.y = (global_ubo.invP * clip_pos).y;

	float xangle = (screen_pos.x - 0.5) * cylindrical_hfov;
	view_dir.x = sin(xangle);
	view_dir.z = cos(xangle);

	view_dir = normalize(view_dir);

	return view_dir * distance;
}

// ---------- Equirectangular transforms ----------
bool equirectangular_forward(vec3 view_pos, out vec2 screen_pos, out float distance)
{
	float lat, lon;
	view_to_latlon(view_pos, lat, lon);
	
	screen_pos.x = (lon / M_PI)*0.5 + 0.5;
	screen_pos.y = (lat / (M_PI/2))*0.5 + 0.5;

	distance = length(view_pos);
	return false;
}

vec3 equirectangular_reverse(vec2 screen_pos, float distance)
{
	//vec3 view_dir = vec3(0, 0, 1);
	vec3 view_dir;
	vec2 angle = (screen_pos * 2 - 1) * vec2(M_PI, M_PI/2);			// Change FOV, change the PI values
	
	latlon_to_view(angle.y, angle.x, view_dir);
	
	return view_dir * distance;
}

// ---------- Mercator transforms ----------
bool mercator_forward(vec3 view_pos, out vec2 screen_pos, out float distance)
{
	vec2 scale = vec2(M_PI, M_PI/2);			// FOV
	float lat, lon;
	view_to_latlon(view_pos, lat, lon);
	screen_pos.x = lon;
	screen_pos.y = log(tan(M_PI*0.25 + lat*0.5));

	screen_pos = screen_pos / scale *0.5 + 0.5;

	distance = length(view_pos);
	return false;
}

vec3 mercator_reverse(vec2 screen_pos, float distance)
{
	vec2 scale = vec2(M_PI, M_PI/2);
	vec3 view_dir;
	vec2 pos = (screen_pos * 2 - 1) * scale;
	//float x = (screen_pos.x * 2 - 1) * M_PI;
	//float y = (screen_pos.y * 2 - 1) * M_PI/2;

	float lon = pos.x;
	float lat = atan(sinh(pos.y));
	
	latlon_to_view(lat, lon, view_dir);
	return view_dir * distance;
}

// ---------- Hammer transforms ----------
	// http://paulbourke.net/geometry/transformationprojection/
vec3 hammer_reverse(vec2 screen_pos, float distance)
{
	vec3 view_dir = vec3(0, 0, 1);
	
	float x = (screen_pos.x * 2 - 1);
	float y = (screen_pos.y * 2 - 1);
	float z = sqrt(1 - x * x / 2 - y * y / 2);
	float longitude = 2 * atan((sqrt(2.0) * z * x)/(2 * z * z - 1));
	float latitude = asin(sqrt(2.0) * z * y);
	
	if (x*x + y*y > 1) {
		//discard;					// Doesn't work
		return vec3(0, 1, 0); 		//TODO set to black, now everyting outside dome is pointing up to ceiling
	}
	
	view_dir.y = -sin(latitude)*view_dir.z;
	view_dir.z = cos(latitude)*view_dir.z;
	
	view_dir.x = sin(longitude)*view_dir.z;
	view_dir.z = cos(longitude)*view_dir.z;
	
	return view_dir * distance;
}

// ---------- Panini transforms ----------
	// http://github.com/shaunlebron/flex-fov
vec2 get_panini(float lat, float lon, float distance)
{
	float d = distance;
	float S = (d+1)/(d+cos(lon));
	float x = S*sin(lon);
	float y = S*tan(lat);
	return vec2(x,y);
}

vec3 get_panini_reverse(vec2 screen_pos, float scale, float distance)
{
	vec3 view_dir;
	screen_pos = (screen_pos * 2 - 1) * scale;
	
	float x = screen_pos.x;
	float y = screen_pos.y/*FIXME *fovy/fovx*/;
	float d = distance;

	float k = x*x/((d+1)*(d+1));
	float dscr = k*k*d*d - (k+1)*(k*d*d-1);
	float clon = (-k*d+sqrt(dscr))/(k+1);
	float S = (d+1)/(d+clon);
	float lon = atan(x,S*clon);
	float lat = atan(y,S);

	view_dir.x = sin(lon)*cos(lat);
	view_dir.y = -sin(lat);
	view_dir.z = cos(lon)*cos(lat);

	return view_dir * distance;
}

// ---------- Ray Direction to Screen Position ----------
bool projection_view_to_screen(vec3 view_pos, out vec2 screen_pos, out float distance, bool previous)
{
	float cylindrical_hfov = previous ? global_ubo.cylindrical_hfov_prev : global_ubo.cylindrical_hfov;
	int projection_type = global_ubo.pt_projection;

	if (projection_type > 3) 		// Replace by direct access to "pt_projection" console variable (not global, have to figure out of to access)
		{ projection_type = 1; }
	
	switch(projection_type)
	{
	case 0:	// Rectilinear
		return rectilinear_forward(view_pos, screen_pos, distance, previous);
	case 1:	// Cylindrical
		return cylindrical_forward(view_pos, screen_pos, distance, previous, cylindrical_hfov);
	case 2:
		return equirectangular_forward(view_pos, screen_pos, distance);
	case 3:
		return mercator_forward(view_pos, screen_pos, distance);
//	case 4:
//		return hammer();
//	case 5:
//		return panini();
	}
}

// ---------- Screen Position to Ray Direction----------
vec3 projection_screen_to_view(vec2 screen_pos, float distance, bool previous)
{
	float cylindrical_hfov = previous ? global_ubo.cylindrical_hfov_prev : global_ubo.cylindrical_hfov;
	int projection_type = global_ubo.pt_projection;

	switch(projection_type)
	{
	case 0:
		return rectlinear_reverse(screen_pos, distance, previous);
	case 1:
		return cylindrical_reverse(screen_pos, distance, previous, cylindrical_hfov);
	case 2:
		return equirectangular_reverse(screen_pos, distance);
	case 3:
		return mercator_reverse(screen_pos, distance);
	case 4:
		return hammer_reverse(screen_pos, distance);
	case 5:
		float scale = get_panini(0, radians(fovx)/2, 1).x;
		return get_panini_reverse(screen_pos, scale, distance);
	}
}
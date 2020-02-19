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

// Conversion from direction vector to angular values in radians and the other way around
void view_to_latlon(vec3 view, out float lat, out float lon)
{
   lon = atan(view.x, view.z);
   lat = atan(-view.y, sqrt(view.x*view.x+view.z*view.z));
}

void latlon_to_view(float lat, float lon, out vec3 view)
{
   view.x = sin(lon)*cos(lat);
   view.y = -sin(lat);
   view.z = cos(lon)*cos(lat);
}

vec2 fov_to_scale()
{
	vec2 scale;
	float aspect = global_ubo.projection_aspect_ratio;
	float fov = global_ubo.pt_projection_fov * M_PI/180.0;		// Convert to radians
	int projection_type = global_ubo.pt_projection;
	
	#define weight_start M_PI/2


	if (fov == 0)		// If FOV is set to 0, use default FOV
	{
		scale = vec2(1.0, 1.0);
		return scale;
	}

	switch(projection_type)
	{
	default:
	case 0:	// Rectilinear
		scale.x = tan(fov/2.0);
		scale.y = tan(fov/2.0) / aspect;
		return scale;

	case 1:	// Cylindrical
		#define progressive 3.0		// Absolute min 2.0 to work up to 360 deg hfov (prefferably more than 2.0)
		scale.x = fov/2.0;
		//scale.y = tan(fov/2.0) / aspect;		// Only works up to FOV = ca 160 deg vertical)
		scale.y = tan(fov/2.0/progressive) * progressive / aspect;
		return scale;

	case 2:	// Equirectangular
		// scale = vec2(M_PI, M_PI/2);	// Full projection
		scale.x = fov/2.0;
		if (fov < weight_start)
		{	scale.y = (fov/2.0) / aspect; }
		else
		{ 
			//scale.y = (fov/2.0) / (cos((fov - M_PI/2)/4*M_PI/8)*(aspect - 2) + 2); 
			scale.y = (fov/2.0) / (  (fov - weight_start)/(2*M_PI-weight_start) * (2 - aspect) + aspect  );
		}
		return scale;

	case 3:
		// vec2 scale = vec2(M_PI, M_PI/2);	// Full projection
		scale.x = fov/2.0;
		if (fov < M_PI/2)
			{ scale.y = (fov/2.0) / aspect; }
		else
		{
			scale.y = (fov/2.0) / (  (fov - weight_start)/(2*M_PI-weight_start) * (2 - aspect) + aspect  );
		}
		return scale;

	case 4:
		return vec2(1,1);

	case 5:
		return vec2(1,1);
	}
}
							   
// ---------- Rectilinear transforms ----------
												  
/*
// Original rectilinear functions (works with build-in "fov" command which assumes 4:3 aspect ratio among other things)
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
*/

bool rectilinear_forward(vec3 view_pos, out vec2 screen_pos, out float distance, bool previous)
{
	vec2 scale = fov_to_scale();
	distance = length(view_pos);
	view_pos = normalize(view_pos);

	float x = view_pos.x;
	float y = -view_pos.y;
	float z = view_pos.z;
	float theta = acos(z);

	if (theta == 0.0)
	{ 	screen_pos = vec2(0.5, 0.5); }
	else
	{
		float r = tan(theta);
		float c = r/sqrt(x*x+y*y);
		screen_pos.x = (x*c/scale.x + 1.0) / 2.0;
		screen_pos.y = (y*c/scale.y + 1.0) / 2.0;
	}

	return screen_pos.y > 0 && screen_pos.y < 1 && screen_pos.x > 0 && screen_pos.x < 1 && view_pos.z > 0;
}

vec3 rectlinear_reverse(vec2 screen_pos, float distance, bool previous)
{
	vec2 scale = fov_to_scale();
	vec3 view_dir;

	float x = (screen_pos.x * 2.0 - 1.0) * scale.x;
	float y = (screen_pos.y * 2.0 - 1.0) * scale.y;

	float r = sqrt(x*x+y*y);
	float theta = atan(r);
	float s = sin(theta);

	view_dir.x = x/r*s;
	view_dir.y = -y/r*s;
	view_dir.z = cos(theta);
	view_dir = normalize(view_dir);

	return view_dir * distance;
}


// ---------- Cylindrical transforms ----------
bool cylindrical_forward(vec3 view_pos, out vec2 screen_pos, out float distance, bool previous)
{
	vec2 scale = fov_to_scale();
	float lat, lon;
	distance = length(view_pos);
	// view_pos = normalize(view_pos);

	view_to_latlon(view_pos, lat, lon);

	float x = lon;
	float y = tan(lat);
	
	screen_pos.x = (x/scale.x + 1.0) / 2.0;
	screen_pos.y = (y/scale.y + 1.0) / 2.0;

	return screen_pos.y > 0 && screen_pos.y < 1 && screen_pos.x > 0 && screen_pos.x < 1;
}

vec3 cylindrical_reverse(vec2 screen_pos, float distance, bool previous)
{
	vec2 scale = fov_to_scale();
	vec3 view_dir;

	float x = (screen_pos.x * 2.0 - 1.0) * scale.x;
	float y = (screen_pos.y * 2.0 - 1.0) * scale.y;

	float lon = x;
	float lat = atan(y);

	latlon_to_view(lat, lon, view_dir);

	return view_dir * distance;
}

// ---------- Equirectangular transforms ----------
bool equirectangular_forward(vec3 view_pos, out vec2 screen_pos, out float distance)
{
	vec2 scale = fov_to_scale();
	float lat, lon;
	distance = length(view_pos);
	// view_pos = normalize(view_pos);

	view_to_latlon(view_pos, lat, lon);
	
	screen_pos.x = (lon / scale.x)*0.5 + 0.5;
	screen_pos.y = (lat / scale.y)*0.5 + 0.5;

	return false;
}

vec3 equirectangular_reverse(vec2 screen_pos, float distance)
{
	vec2 scale = fov_to_scale();
	vec3 view_dir;
	vec2 angle = (screen_pos * 2 - 1) * scale;			// Change FOV, change the PI values
	
	latlon_to_view(angle.y, angle.x, view_dir);
	
	return view_dir * distance;
}

// ---------- Mercator transforms ----------
bool mercator_forward(vec3 view_pos, out vec2 screen_pos, out float distance)
{
	vec2 scale = fov_to_scale();
	float lat, lon;
	distance = length(view_pos);
	view_pos = normalize(view_pos);

	view_to_latlon(view_pos, lat, lon);
	screen_pos.x = lon;
	screen_pos.y = log(tan(M_PI*0.25 + lat*0.5));

	screen_pos = screen_pos / scale *0.5 + 0.5;

	return false;
}

vec3 mercator_reverse(vec2 screen_pos, float distance)
{
	vec2 scale = fov_to_scale();
	vec3 view_dir;
	vec2 pos = (screen_pos * 2 - 1) * scale;

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
		return vec3(0, 0, 0); 		//TODO set to black, now everyting outside reflects sky color
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
	float y = screen_pos.y;
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

// ---------- Temporal Filtering: Ray Direction to Screen Position ----------
bool projection_view_to_screen(vec3 view_pos, out vec2 screen_pos, out float distance, bool previous)
{
	int projection_type = global_ubo.pt_projection;
	float aspect = global_ubo.projection_aspect_ratio;
	float fov = global_ubo.pt_projection_fov * M_PI/180;		// Convert to radians

	if (projection_type > 3) 		// Replace by direct access to "pt_projection" console variable (not global, have to figure out of to access)
		{ projection_type = 1; }
	
	switch(projection_type)
	{
	case 0:	// Rectilinear
		return rectilinear_forward(view_pos, screen_pos, distance, previous);
	case 1:	// Cylindrical
		return cylindrical_forward(view_pos, screen_pos, distance, previous);
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

// ---------- Rendering: Screen Position to Ray Direction----------
vec3 projection_screen_to_view(vec2 screen_pos, float distance, bool previous)
{
	int projection_type = global_ubo.pt_projection;
	float aspect = global_ubo.projection_aspect_ratio;
	float fov = global_ubo.pt_projection_fov * M_PI/180;		// Convert to radians

	switch(projection_type)
	{
	default:
	case 0:
		return rectlinear_reverse(screen_pos, distance, previous);
	case 1:
		return cylindrical_reverse(screen_pos, distance, previous);
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
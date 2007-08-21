/*
Copyright (C) 2003-2006 Andrey Nazarov

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

//
// field.c
//

#include "config.h"
#include "q_shared.h"
#include "com_public.h"
#include "key_public.h"
#include "ref_public.h"
#include "q_field.h"

/*
================
IF_Init
================
*/
void IF_Init( inputField_t *field, int visibleChars, int maxChars ) {
	memset( field, 0, sizeof( *field ) );
	
	clamp( maxChars, 1, sizeof( field->text ) );
	clamp( visibleChars, 1, maxChars );

	field->maxChars = maxChars;
	field->visibleChars = visibleChars;
}

/*
================
IF_InitText
================
*/
void IF_InitText( inputField_t *field, int visibleChars, int maxChars,
        const char *text )
{
	memset( field, 0, sizeof( *field ) );
	
	clamp( maxChars, 1, sizeof( field->text ) );
	clamp( visibleChars, 1, maxChars );

	field->maxChars = maxChars;
	field->visibleChars = visibleChars;

	Q_strncpyz( field->text, text, sizeof( field->text ) );
	field->cursorPos = strlen( field->text );
}

/*
================
IF_Clear
================
*/
void IF_Clear( inputField_t *field ) {
	memset( field->text, 0, sizeof( field->text ) );
	field->cursorPos = 0;
	field->selectStart = 0;
	field->selectEnd = 0;
}

/*
================
IF_Replace
================
*/
void IF_Replace( inputField_t *field, const char *text ) {
	Q_strncpyz( field->text, text, sizeof( field->text ) );
	field->cursorPos = strlen( field->text );
	field->selectStart = 0;
	field->selectEnd = 0;
}

#ifndef DEDICATED_ONLY


/*
================
IF_DeleteSelection
================
*/
void IF_DeleteSelection( inputField_t *field ) {
	if( field->selectStart < field->selectEnd ) {
		memmove( field->text + field->selectStart,
                field->text + field->selectEnd,
                    sizeof( field->text ) - field->selectStart );
		field->selectEnd = field->selectStart;
	}
}

/*
================
IF_KeyEvent
================
*/
qboolean IF_KeyEvent( inputField_t *field, int key ) {
	if( key == K_DEL ) {
		if( field->selectStart < field->selectEnd ) {
			IF_DeleteSelection( field );
			return qtrue;
		}
		if( field->text[field->cursorPos] ) {
			memmove( field->text + field->cursorPos,
                field->text + field->cursorPos + 1,
                sizeof( field->text ) - field->cursorPos );
		}
		return qtrue;
	}
	
	if( key == K_BACKSPACE || ( key == 'h' && keys.IsDown( K_CTRL ) ) ) {
		if( field->selectStart < field->selectEnd ) {
			IF_DeleteSelection( field );
			return qtrue;
		}
		if( field->cursorPos > 0 ) {
			memmove( field->text + field->cursorPos - 1,
                field->text + field->cursorPos,
                sizeof( field->text ) - field->cursorPos );
			field->cursorPos--;
		}
		return qtrue;
	}

    if( key == 'w' && keys.IsDown( K_CTRL ) ) {
        int oldpos = field->cursorPos;

        // kill trailing whitespace
		while( field->cursorPos > 0 ) {
            if( field->text[field->cursorPos] > ' ' ) {
                break;
            }
			field->cursorPos--;
        }

        // kill this word
		while( field->cursorPos > 0 ) {
            if( field->text[ field->cursorPos - 1 ] <= ' ' ) {
                break;
            }
			field->cursorPos--;
        }
		memmove( field->text + field->cursorPos ,
                field->text + oldpos,
                sizeof( field->text ) - oldpos );
    }

    if( key == 'k' && keys.IsDown( K_CTRL ) ) {
        IF_Clear( field );
        return qtrue;
    }

	if( ( ( key == 'V' || key == 'v' ) && keys.IsDown( K_CTRL ) ) ||
		( key == K_INS && keys.IsDown( K_SHIFT ) ) ||
        key == K_MOUSE3 )
	{
		char *cbd, *s;
		
		if( ( cbd = sys.GetClipboardData() ) != NULL ) {
			s = cbd;
			while( *s ) {
                switch( *s ) {
                case '\n':
                    if( s[1] ) {
    					IF_CharEvent( field, ';' );
                    }
                    break;
                case '\r':
                case '\t':
    			    IF_CharEvent( field, ' ' );
                    break;
                default:
					IF_CharEvent( field, *s );
                    break;
				}
				s++;
			}
			com.Free( cbd );
		}

		return qtrue;
	}

    if( key == 'c' && keys.IsDown( K_CTRL ) ) {
        sys.SetClipboardData( field->text );
        return qtrue;
    }
	
	if( key == K_LEFTARROW ) {
		if( field->cursorPos > 0 ) {
			field->cursorPos--;
		}
		return qtrue;
	}

	if( key == K_RIGHTARROW ) {
		if( field->text[field->cursorPos] ) {
			field->cursorPos++;
		}
		return qtrue;
	}

	if( key == K_HOME ) {
		field->cursorPos = 0;
		return qtrue;
	}

	if( key == K_END ) {
		field->cursorPos = strlen( field->text );
		return qtrue;
	}

	if( key == K_INS ) {
		keys.SetOverstrikeMode( keys.GetOverstrikeMode() ^ 1 );
		return qtrue;
	}

	return qfalse;
}

/*
================
IF_CharEvent
================
*/
qboolean IF_CharEvent( inputField_t *field, int key ) {
	if( key < 32 || key > 127 ) {
		return qfalse;	// non printable
	}

	if( field->cursorPos >= field->maxChars - 1 ) {
		// buffer limit was reached, just replace the last character
		field->text[field->cursorPos] = key;
		return qtrue;
	}

	if( field->selectStart < field->selectEnd ) {
		IF_DeleteSelection( field );
		field->cursorPos = field->selectStart;
		field->text[field->cursorPos] = key;
		return qtrue;
		
	}

	if( keys.GetOverstrikeMode() ) {
		// replace the character at cursor and advance
		field->text[field->cursorPos] = key;
		field->cursorPos++;
		return qtrue;
	}

	// insert new character at cursor position
	memmove( field->text + field->cursorPos + 1,
            field->text + field->cursorPos,
                sizeof( field->text ) - field->cursorPos - 1 );
	field->text[field->cursorPos] = key;
	field->cursorPos++;

	return qtrue;
}

/*
================
IF_Draw

The input line scrolls horizontally if typing goes beyond the right edge
================
*/
void IF_Draw( inputField_t *field, int x, int y, uint32 flags,
        qhandle_t hFont )
{
	char *text;
	int cursorPos, offset;
	int start, width;
	int cw, ch;
	color_t	selectColor;

	if( field->cursorPos > sizeof( field->text ) - 1 ) {
		Com_Error( ERR_FATAL, "IF_Draw: bad field->cursorPos" );
	}

	text = field->text;
	cursorPos = field->cursorPos;
	offset = 0;

	/* scroll horizontally */
	if( cursorPos > field->visibleChars - 1 ) {
		cursorPos = field->visibleChars - 1;
		offset = field->cursorPos - cursorPos;
	}

	if( !( flags & UI_DRAWCURSOR ) ) {
        /* just draw text and return */
        ref.DrawString( x, y, flags,
                field->visibleChars, text + offset, hFont );
        return;
    }
    
	ref.DrawGetFontSize( &cw, &ch, hFont );

    /* draw selection background */
    if( field->selectStart < field->selectEnd ) {
        start = field->selectStart - offset;
        width = field->selectEnd - field->selectStart;

        if( width > ( field->visibleChars - 1 ) - start ) {
            width = ( field->visibleChars - 1 ) - start;
        }

        Vector4Set( selectColor, 156, 90, 84, 255 );
        ref.DrawFillEx( x + start * cw, y, width * cw, ch, selectColor );
    }

	/* draw text */
	ref.DrawString( x, y, flags,
            field->visibleChars, text + offset, hFont );

    /* draw blinking cursor */
    if( ( sys.Milliseconds() >> 8 ) & 1 ) {
        int c = keys.GetOverstrikeMode() ? 11 : '_';
        ref.DrawChar( x + cursorPos * cw, y, flags, c, hFont );
    }

}

#endif
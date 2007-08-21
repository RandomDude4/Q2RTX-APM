/*
Copyright (C) 1997-2001 Id Software, Inc.

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

#include "com_local.h"
#include "net_chan.h"

/*

packet header
-------------
31	sequence
1	does this message contain a reliable payload
31	acknowledge sequence
1	acknowledge receipt of even/odd message
16	qport

The remote connection never knows if it missed a reliable message, the
local side detects that it has been dropped by seeing a sequence acknowledge
higher thatn the last reliable sequence, but without the correct evon/odd
bit for the reliable set.

If the sender notices that a reliable message has been dropped, it will be
retransmitted.  It will not be retransmitted again until a message after
the retransmit has been acknowledged and the reliable still failed to get theref.

if the sequence number is -1, the packet should be handled without a netcon

The reliable message can be added to at any time by doing
MSG_Write* (&netchan->message, <data>).

If the message buffer is overflowed, either by a single message, or by
multiple frames worth piling up while the last reliable transmit goes
unacknowledged, the netchan signals a fatal error.

Reliable messages are always placed first in a packet, then the unreliable
message is included if there is sufficient room.

To the receiver, there is no distinction between the reliable and unreliable
parts of the message, they are just processed out as a single larger message.

Illogical packet sequence numbers cause the packet to be dropped, but do
not kill the connection.  This, combined with the tight window of valid
reliable acknowledgement numbers provides protection against malicious
address spoofing.


The qport field is a workaround for bad address translating routers that
sometimes remap the client's source port on a packet during gameplay.

If the base part of the net address matches and the qport matches, then the
channel matches even if the IP port differs.  The IP port should be updated
to the new value before sending out any replies.


If there is no information that needs to be transfered on a given frame,
such as during the connection stage while waiting for the client to load,
then a packet only needs to be delivered if there is something in the
unacknowledged reliable
*/

static cvar_t		*showpackets;
static cvar_t		*showdrop;

cvar_t		*net_qport;
cvar_t      *net_maxmsglen;
cvar_t      *net_chantype;

netadr_t	net_from;

/*
===============
Netchan_Init

===============
*/
void Netchan_Init( void ) {
	int		port;

	// pick a port value that should be nice and random
	port = Sys_Milliseconds() & 0xffff;

	showpackets = Cvar_Get( "showpackets", "0", 0 );
	showdrop = Cvar_Get( "showdrop", "0", 0 );
	net_qport = Cvar_Get( "qport", va( "%i", port ), 0 );
    net_maxmsglen = Cvar_Get( "net_maxmsglen", "1390", 0 );
    net_chantype = Cvar_Get( "net_chantype", "1", 0 );
}

/*
===============
Netchan_OutOfBand

Sends an out-of-band datagram
================
*/
neterr_t Netchan_OutOfBand( netsrc_t sock, const netadr_t *address,
        uint32 length, const byte *data )
{
	sizebuf_t	send;
	byte		send_data[MAX_PACKETLEN_DEFAULT];
	neterr_t	ret;

	SZ_Init( &send, send_data, sizeof( send_data ) );

// write the packet header
	SZ_WriteLong( &send, -1 );	// -1 sequence means out of band
	SZ_Write( &send, data, length );

// send the datagram
	ret = NET_SendPacket( sock, address, send.cursize, send.data );

	return ret;

}

/*
===============
Netchan_OutOfBandPrint

Sends a text message in an out-of-band datagram
================
*/
neterr_t Netchan_OutOfBandPrint( netsrc_t sock, const netadr_t *address,
        const char *format, ... )
{
	va_list		argptr;
	byte		send_data[MAX_PACKETLEN_DEFAULT];
	int			length;
	neterr_t	ret;

/* write the packet header */
	*( uint32 * )send_data = -1;	/* -1 sequence means out of band */
	
	va_start( argptr, format );
	length = Q_vsnprintf( ( char * )send_data + 4, sizeof( send_data ) - 4,
            format, argptr );
	va_end( argptr );

/* send the datagram */
	ret = NET_SendPacket( sock, address, length + 4, send_data );

	return ret;
}

// ============================================================================

static int NetchanOld_TransmitNextFragment( netchan_t *netchan ) {
	Com_Error( ERR_FATAL, "NetchanOld_TransmitNextFragment: not implemented" );
	return 0;
}

/*
===============
NetchanOld_Transmit

tries to send an unreliable message to a connection, and handles the
transmition / retransmition of the reliable messages.

A 0 length will still generate a packet and deal with the reliable messages.
================
*/
static int NetchanOld_Transmit( netchan_t *netchan, int length,
        const byte *data )
{
	netchan_old_t *chan = ( netchan_old_t * )netchan;
	sizebuf_t	send;
	byte		send_buf[MAX_PACKETLEN];
	qboolean	send_reliable;
	uint32	    w1, w2;
	neterr_t	ret;

// check for message overflow
	if( netchan->message.overflowed ) {
		netchan->fatal_error = qtrue;
		Com_WPrintf( "%s: outgoing message overflow\n",
                NET_AdrToString( &netchan->remote_address ) );
		return -1;
	}

	send_reliable = qfalse;

	// if the remote side dropped the last reliable message, resend it
	if( netchan->incoming_acknowledged > chan->last_reliable_sequence &&
		chan->incoming_reliable_acknowledged != chan->reliable_sequence )
	{
		send_reliable = qtrue;
	}

// if the reliable transmit buffer is empty, copy the current message out
	if( !netchan->reliable_length && netchan->message.cursize ) {
		send_reliable = qtrue;
		memcpy( chan->reliable_buf, chan->message_buf,
                netchan->message.cursize );
		netchan->reliable_length = netchan->message.cursize;
		netchan->message.cursize = 0;
		chan->reliable_sequence ^= 1;
	}

// write the packet header
	w1 = ( netchan->outgoing_sequence & ~( 1 << 31 ) ) |
        ( send_reliable << 31 );
	w2 = ( netchan->incoming_sequence & ~( 1 << 31 ) ) |
        ( chan->incoming_reliable_sequence << 31 );

	netchan->outgoing_sequence++;
	netchan->reliable_ack_pending = qfalse;
	netchan->last_sent = com_localTime;

	SZ_Init( &send, send_buf, sizeof( send_buf ) );

	SZ_WriteLong( &send, w1 );
	SZ_WriteLong( &send, w2 );

	// send the qport if we are a client
	if( netchan->sock == NS_CLIENT ) {
		if( netchan->protocol < PROTOCOL_VERSION_R1Q2 ) {
			SZ_WriteShort( &send, netchan->qport );
		} else if( netchan->qport ) {
			SZ_WriteByte( &send, netchan->qport );
		}
	}

// copy the reliable message to the packet first
	if( send_reliable ) {
		SZ_Write( &send, chan->reliable_buf, netchan->reliable_length );
		chan->last_reliable_sequence = netchan->outgoing_sequence;
	}
	
// add the unreliable part if space is available
	if( send.maxsize - send.cursize >= length )
		SZ_Write( &send, data, length );
	else
		Com_WPrintf( "%s: dumped unreliable\n",
            NET_AdrToString( &netchan->remote_address ) );

	if( showpackets->integer ) {
        Com_Printf( "send %4i : s=%i ack=%i rack=%i",
            send.cursize,
            netchan->outgoing_sequence - 1,
            netchan->incoming_sequence,
            chan->incoming_reliable_sequence );
		if( send_reliable ) {
			Com_Printf( " reliable=%i", chan->reliable_sequence );
		}
        Com_Printf( "\n" );
	}

	// send the datagram
	ret = NET_SendPacket( netchan->sock, &netchan->remote_address,
        send.cursize, send.data );
	if( ret == NET_ERROR ) {
		return -1;
	}

	return send.cursize;

}

/*
=================
NetchanOld_Process

called when the current net_message is from remote_address
modifies net_message so that it points to the packet payload
=================
*/
static qboolean NetchanOld_Process( netchan_t *netchan ) {
	netchan_old_t *chan = ( netchan_old_t * )netchan;
	uint32  	sequence, sequence_ack;
	uint32  	reliable_ack, reliable_message;
	int			qport;

// get sequence numbers		
	MSG_BeginReading();
	sequence = MSG_ReadLong();
	sequence_ack = MSG_ReadLong();

	// read the qport if we are a server
	if( netchan->sock == NS_SERVER ) {
		if( netchan->protocol < PROTOCOL_VERSION_R1Q2 ) {
			qport = MSG_ReadShort();
		} else if( netchan->qport ) {
			qport = MSG_ReadByte();
		}
	}

	reliable_message = sequence >> 31;
	reliable_ack = sequence_ack >> 31;

	sequence &= ~( 1 << 31 );
	sequence_ack &= ~( 1 << 31 );	

	if( showpackets->integer ) {
        Com_Printf( "recv %4i : s=%i ack=%i rack=%i",
            msg_read.cursize,
            sequence,
            sequence_ack,
            reliable_ack );
		if( reliable_message ) {
			Com_Printf( " reliable=%i", chan->incoming_reliable_sequence ^ 1 );
        }
        Com_Printf( "\n" );
	}

//
// discard stale or duplicated packets
//
	if( sequence <= netchan->incoming_sequence ) {
		if( showdrop->integer )
			Com_Printf( "%s: out of order packet %i at %i\n",
				NET_AdrToString( &netchan->remote_address ),
				sequence,
				netchan->incoming_sequence );
		return qfalse;
	}

//
// dropped packets don't keep the message from being used
//
	netchan->dropped = sequence - ( netchan->incoming_sequence + 1 );
	if( netchan->dropped > 0 ) {
		if( showdrop->integer )
			Com_Printf( "%s: dropped %i packets at %i\n",
			NET_AdrToString( &netchan->remote_address ),
			netchan->dropped,
			sequence );
	}

//
// if the current outgoing reliable message has been acknowledged
// clear the buffer to make way for the next
//
	chan->incoming_reliable_acknowledged = reliable_ack;
	if( reliable_ack == chan->reliable_sequence )
		netchan->reliable_length = 0;	// it has been received
	
//
// if this message contains a reliable message, bump incoming_reliable_sequence 
//
	netchan->incoming_sequence = sequence;
	netchan->incoming_acknowledged = sequence_ack;
	if( reliable_message ) {
		netchan->reliable_ack_pending = qtrue;
		chan->incoming_reliable_sequence ^= 1;
	}

//
// the message can now be read from the current message pointer
//
	netchan->last_received = com_localTime;

	return qtrue;
}

/*
===============
NetchanOld_ShouldUpdate
================
*/
static qboolean NetchanOld_ShouldUpdate( netchan_t *netchan ) {
	if( netchan->message.cursize || netchan->reliable_ack_pending ||
            com_localTime - netchan->last_sent > 1000 )
    {
		return qtrue;
	}

	return qfalse;
}

/*
==============
NetchanOld_Setup

called to open a channel to a remote system
==============
*/
static netchan_t *NetchanOld_Setup( netsrc_t sock, const netadr_t *adr,
        int qport, int maxpacketlen )
{
	netchan_old_t *chan;
	netchan_t *netchan;

	Z_Reserve( sizeof( *chan ) + maxpacketlen * 2 );

    chan = Z_ReservedAlloc( sizeof( *chan ) );
    memset( chan, 0, sizeof( *chan ) );
	netchan = ( netchan_t * )chan;
	netchan->sock = sock;
	netchan->remote_address = *adr;
	netchan->qport = qport;
    netchan->maxpacketlen = maxpacketlen;
	netchan->last_received = com_localTime;
	netchan->incoming_sequence = 0;
	netchan->outgoing_sequence = 1;

	netchan->Process = NetchanOld_Process;
	netchan->Transmit = NetchanOld_Transmit;
	netchan->TransmitNextFragment = NetchanOld_TransmitNextFragment;
	netchan->ShouldUpdate = NetchanOld_ShouldUpdate;

    chan->message_buf = Z_ReservedAlloc( maxpacketlen );
	SZ_Init( &netchan->message, chan->message_buf, maxpacketlen );

    chan->reliable_buf = Z_ReservedAlloc( maxpacketlen );

	return netchan;
}

// ============================================================================


/*
===============
NetchanNew_TransmitNextFragment
================
*/
static int NetchanNew_TransmitNextFragment( netchan_t *netchan ) {
	netchan_new_t *chan = ( netchan_new_t * )netchan;
	sizebuf_t	send;
	byte		send_buf[MAX_PACKETLEN];
	qboolean	send_reliable;
	uint32	    w1, w2;
	uint16	    offset;
	int			fragment_length;
	qboolean	more_fragments;
	neterr_t	ret;

	send_reliable = netchan->reliable_length ? qtrue : qfalse;

	/* write the packet header */
	w1 = ( netchan->outgoing_sequence & 0x3FFFFFFF ) | ( 1 << 30 ) |
        ( send_reliable << 31 );
	w2 = ( netchan->incoming_sequence & 0x3FFFFFFF ) | ( 0 << 30 ) |
        ( chan->incoming_reliable_sequence << 31 );

	SZ_Init( &send, send_buf, sizeof( send_buf ) );
	send.allowoverflow = qfalse;

	SZ_WriteLong( &send, w1 );
	SZ_WriteLong( &send, w2 );

	/* send the qport if we are a client */
	if( netchan->sock == NS_CLIENT && netchan->qport ) {
		SZ_WriteByte( &send, netchan->qport );
	}

	fragment_length = chan->fragment_out.cursize - chan->fragment_out.readcount;
	if( fragment_length > netchan->maxpacketlen ) {
		fragment_length = netchan->maxpacketlen;
	}
	
	more_fragments = qtrue;
	if( chan->fragment_out.readcount + fragment_length ==
        chan->fragment_out.cursize )
    {
		more_fragments = qfalse;
	}

	/* write fragment offset */
	offset = ( chan->fragment_out.readcount & 0x7FFF ) |
        ( more_fragments << 15 );
	SZ_WriteShort( &send, offset );

	/* write fragment contents */
	SZ_Write( &send, chan->fragment_out.data + chan->fragment_out.readcount,
            fragment_length );
	
	if( showpackets->integer ) {
		Com_Printf( "send %4i : s=%i ack=%i rack=%i "
                    "fragment_offset=%i more_fragments=%i",
			send.cursize,
			netchan->outgoing_sequence,
			netchan->incoming_sequence,
			chan->incoming_reliable_sequence,
			chan->fragment_out.readcount,
			more_fragments );
		if( send_reliable ) {
			Com_Printf( " reliable=%i ", chan->reliable_sequence );
		}
		Com_Printf( "\n" );
	}

	chan->fragment_out.readcount += fragment_length;
	netchan->fragment_pending = more_fragments;

	/* if the message has been sent completely, clear the fragment buffer */
	if( !netchan->fragment_pending ) {
		netchan->outgoing_sequence++;
		netchan->last_sent = com_localTime;
		SZ_Clear( &chan->fragment_out );
	}

	// send the datagram
	ret = NET_SendPacket( netchan->sock, &netchan->remote_address,
            send.cursize, send.data );
	if( ret == NET_ERROR ) {
		return -1;
	}

	return send.cursize;
}

/*
===============
NetchanNew_Transmit
================
*/
static int NetchanNew_Transmit( netchan_t *netchan, int length,
        const byte *data )
{
	netchan_new_t *chan = ( netchan_new_t * )netchan;
	sizebuf_t	send;
	byte		send_buf[MAX_PACKETLEN];
	qboolean	send_reliable;
	uint32		w1, w2;
	neterr_t	ret;

// check for message overflow
	if( netchan->message.overflowed ) {
		netchan->fatal_error = qtrue;
		Com_WPrintf( "%s: outgoing message overflow\n",
                NET_AdrToString( &netchan->remote_address ) );
		return -1;
	}

	if( netchan->fragment_pending ) {
		return NetchanNew_TransmitNextFragment( netchan );
	}

	send_reliable = qfalse;

// if the remote side dropped the last reliable message, resend it
	if( netchan->incoming_acknowledged > chan->last_reliable_sequence &&
		chan->incoming_reliable_acknowledged != chan->reliable_sequence )
	{
		send_reliable = qtrue;
	}

// if the reliable transmit buffer is empty, copy the current message out
	if( !netchan->reliable_length && netchan->message.cursize ) {
		send_reliable = qtrue;
		memcpy( chan->reliable_buf, chan->message_buf,
                netchan->message.cursize );
		netchan->reliable_length = netchan->message.cursize;
		netchan->message.cursize = 0;
		chan->reliable_sequence ^= 1;
	}

	if( length > netchan->maxpacketlen || ( send_reliable &&
        ( netchan->reliable_length + length > netchan->maxpacketlen ) ) )
    {
		if( send_reliable ) {
			chan->last_reliable_sequence = netchan->outgoing_sequence;
			SZ_Write( &chan->fragment_out, chan->reliable_buf,
                    netchan->reliable_length );
		}
		SZ_Write( &chan->fragment_out, data, length );
		return NetchanNew_TransmitNextFragment( netchan );
	}

// write the packet header
	w1 = ( netchan->outgoing_sequence & 0x3FFFFFFF ) | ( send_reliable << 31 );
	w2 = ( netchan->incoming_sequence & 0x3FFFFFFF ) |
        ( chan->incoming_reliable_sequence << 31 );

	netchan->outgoing_sequence++;
	netchan->reliable_ack_pending = qfalse;
	netchan->last_sent = com_localTime;

	SZ_Init( &send, send_buf, netchan->maxpacketlen );
	send.allowoverflow = qfalse;

	SZ_WriteLong( &send, w1 );
	SZ_WriteLong( &send, w2 );

	// send the qport if we are a client
	if( netchan->sock == NS_CLIENT && netchan->qport ) {
		SZ_WriteByte( &send, netchan->qport );
	}
	
	// copy the reliable message to the packet first
	if( send_reliable ) {
		chan->last_reliable_sequence = netchan->outgoing_sequence;
		SZ_Write( &send, chan->reliable_buf, netchan->reliable_length );
	}
	
	// add the unreliable part
	SZ_Write( &send, data, length );

	if( showpackets->integer ) {
        Com_Printf( "send %4i : s=%i ack=%i rack=%i",
            send.cursize,
            netchan->outgoing_sequence - 1,
            netchan->incoming_sequence,
            chan->incoming_reliable_sequence );
		if( send_reliable ) {
			Com_Printf( " reliable=%i", chan->reliable_sequence );
        }
        Com_Printf( "\n" );
	}

	// send the datagram
	ret = NET_SendPacket( netchan->sock, &netchan->remote_address,
            send.cursize, send.data );
	if( ret == NET_ERROR ) {
		return -1;
	}

	return send.cursize;

}

/*
=================
NetchanNew_Process
=================
*/
static qboolean NetchanNew_Process( netchan_t *netchan ) {
	netchan_new_t *chan = ( netchan_new_t * )netchan;
	uint32      sequence, sequence_ack, reliable_ack;
    qboolean    reliable_message, fragmented_message, more_fragments;
	uint16		fragment_offset;
	int			length;

// get sequence numbers		
	MSG_BeginReading();
	sequence = MSG_ReadLong();
	sequence_ack = MSG_ReadLong();

	// read the qport if we are a server
	if( netchan->sock == NS_SERVER && netchan->qport ) {
		MSG_ReadByte();
	}

	reliable_message = sequence >> 31;
	reliable_ack = sequence_ack >> 31;
	fragmented_message = ( sequence >> 30 ) & 1;

	sequence &= 0x3FFFFFFF;
	sequence_ack &= 0x3FFFFFFF;

    fragment_offset = 0;
    more_fragments = qfalse;
	if( fragmented_message ) {
		fragment_offset = MSG_ReadShort();
        more_fragments = fragment_offset >> 15;
		fragment_offset &= 0x7FFF;
	}

	if( showpackets->integer ) {
		Com_Printf( "recv %4i : s=%i ack=%i rack=%i", msg_read.cursize,
                sequence, sequence_ack, reliable_ack );
		if( fragmented_message ) {
			Com_Printf( " fragment_offset=%i more_fragments=%i",
                    fragment_offset, more_fragments );
		}
		if( reliable_message ) {
			Com_Printf( " reliable=%i", chan->incoming_reliable_sequence ^ 1 );
		}
		Com_Printf( "\n" );
	}

//
// discard stale or duplicated packets
//
	if( sequence <= netchan->incoming_sequence ) {
		if( showdrop->integer || showpackets->integer ) {
			Com_Printf( "%s: out of order packet %i at %i\n",
                NET_AdrToString( &netchan->remote_address ),
                    sequence, netchan->incoming_sequence );
		}
		return qfalse;
	}

//
// dropped packets don't keep the message from being used
//
	netchan->dropped = sequence - ( netchan->incoming_sequence + 1 );
	if( netchan->dropped > 0 ) {
		if( showdrop->integer || showpackets->integer ) {
			Com_Printf( "%s: dropped %i packets at %i\n",
                NET_AdrToString( &netchan->remote_address ),
                    netchan->dropped, sequence );
		}
	}

//
// if the current outgoing reliable message has been acknowledged
// clear the buffer to make way for the next
//
	chan->incoming_reliable_acknowledged = reliable_ack;
	if( reliable_ack == chan->reliable_sequence ) {
		netchan->reliable_length = 0;	/* it has been received */
	}


//
// parse fragment header, if any
//
	if( fragmented_message ) {
		if( chan->fragment_sequence != sequence ) {
			/* start new receive sequence */
			chan->fragment_sequence = sequence;
			SZ_Clear( &chan->fragment_in );
		}

		if( fragment_offset < chan->fragment_in.cursize ) {
			if( showdrop->integer || showpackets->integer ) {
				Com_Printf( "%s: out of order fragment at %i\n",
                    NET_AdrToString( &netchan->remote_address ),
                        sequence );
			}
			return qfalse;
		}

		if( fragment_offset > chan->fragment_in.cursize ) {
			if( showdrop->integer || showpackets->integer ) {
				Com_Printf( "%s: dropped fragment(s) at %i\n",
                    NET_AdrToString( &netchan->remote_address ),
                        sequence );
			}
			return qfalse;
		}

		length = msg_read.cursize - msg_read.readcount;
		if( chan->fragment_in.cursize + length > chan->fragment_in.maxsize ) {
			if( showdrop->integer || showpackets->integer ) {
				Com_Printf( "%s: oversize fragment at %i\n",
                    NET_AdrToString( &netchan->remote_address ),
                        sequence );
			}
			return qfalse;
		}

		SZ_Write( &chan->fragment_in, msg_read.data + msg_read.readcount,
                length );
		if( more_fragments ) {
			return qfalse;
		}

		/* message has been sucessfully assembled */
		SZ_Clear( &msg_read );
		SZ_Write( &msg_read, chan->fragment_in.data,
                chan->fragment_in.cursize );
		SZ_Clear( &chan->fragment_in );
	}
	
	netchan->incoming_sequence = sequence;
	netchan->incoming_acknowledged = sequence_ack;
	
//
// if this message contains a reliable message, bump incoming_reliable_sequence 
//
	if( reliable_message ) {
		netchan->reliable_ack_pending = qtrue;
		chan->incoming_reliable_sequence ^= 1;
	}

//
// the message can now be read from the current message pointer
//
	netchan->last_received = com_localTime;

	return qtrue;
}

/*
==============
NetchanNew_ShouldUpdate
==============
*/
static qboolean NetchanNew_ShouldUpdate( netchan_t *netchan ) {
	netchan_new_t *chan = ( netchan_new_t * )netchan;

	if( netchan->message.cursize ||
        netchan->reliable_ack_pending ||
		chan->fragment_out.cursize ||
        com_localTime - netchan->last_sent > 1000 )
	{
		return qtrue;
	}

	return qfalse;
}

/*
==============
NetchanNew_Setup
==============
*/
static netchan_t *NetchanNew_Setup( netsrc_t sock, const netadr_t *adr,
        int qport, int maxpacketlen )
{
	netchan_new_t *chan;
	netchan_t *netchan;

	chan = Z_Mallocz( sizeof( *chan ) );
	netchan = ( netchan_t * )chan;
	netchan->sock = sock;
	netchan->remote_address = *adr;
	netchan->qport = qport;
    netchan->maxpacketlen = maxpacketlen;
	netchan->last_received = com_localTime;
	netchan->incoming_sequence = 0;
	netchan->outgoing_sequence = 1;

	netchan->Process = NetchanNew_Process;
	netchan->Transmit = NetchanNew_Transmit;
	netchan->TransmitNextFragment = NetchanNew_TransmitNextFragment;
	netchan->ShouldUpdate = NetchanNew_ShouldUpdate;

	SZ_Init( &netchan->message, chan->message_buf,
            sizeof( chan->message_buf ) );
	SZ_Init( &chan->fragment_in, chan->fragment_in_buf,
            sizeof( chan->fragment_in_buf ) );
	SZ_Init( &chan->fragment_out, chan->fragment_out_buf,
            sizeof( chan->fragment_out_buf ) );

	return netchan;
}

/*
==============
Netchan_Setup
==============
*/
netchan_t *Netchan_Setup( netsrc_t sock, netchan_type_t type,
        const netadr_t *adr, int qport, int maxpacketlen, int protocol )
{
    netchan_t *netchan;

    clamp( maxpacketlen, 256, MAX_PACKETLEN_WRITABLE );

	switch( type ) {
	case NETCHAN_OLD:
		netchan = NetchanOld_Setup( sock, adr, qport, maxpacketlen );
        break;
	case NETCHAN_NEW:
		netchan = NetchanNew_Setup( sock, adr, qport, maxpacketlen );
        break;
	default:
		Com_Error( ERR_FATAL, "Netchan_Setup: bad type" );
        netchan = NULL;
	}

    netchan->protocol = protocol;
    netchan->type = type;

	return netchan;

}

/*
==============
Netchan_Close
==============
*/
void Netchan_Close( netchan_t *netchan ) {
	Z_Free( netchan );
}

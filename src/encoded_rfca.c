/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include "encoded_rfca.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// We use this value to temporary mask or select values in the rfca,
// under the assumption that a base this large is unfeasible and will never occur.
const rfca_node_t MASKED_VALUE = (1 << 31);

pattern_list_t*
standardCodeTable( int base ) {
    pattern_list_t* ct;
    ct = (pattern_list_t*)malloc( sizeof( pattern_list_t ) );
    INIT_LIST_HEAD( &(ct->list) );

    for( int i =0; i < base; i++ ) {
        pattern_list_t* tmp;
        tmp = (pattern_list_t*)malloc( sizeof( pattern_list_t ) );
        INIT_LIST_HEAD( &(tmp->list) );

        tmp->pattern = pattern_createSingle( i );
        list_add_tail( &(tmp->list), &(ct->list ) );
    }
    return ct;
}

double
computeCodeTableBits( encoded_rfca_t* v ) {
    double bits =0.0;
    struct list_head* pos;
    list_for_each( pos, &(v->codeTable->list) ) {
        pattern_list_t* tmp = list_entry( pos, pattern_list_t, list );
        bits += tmp->pattern->codeLength + v->stdBitsPerSingleton * tmp->pattern->size;
    }
    return bits;
}

double
computeEncodedBits( encoded_rfca_t* v ) {
    double bits =0.0;
    struct list_head* pos;
    list_for_each( pos, &(v->encoded->list) ) {
        region_list_t* tmp = list_entry( pos, region_list_t, list );
        region_t* region = tmp->region;
        bits += region->pattern->codeLength;
    }
    return bits;
}

void
maskRegion( rfca_t* r, pattern_t* p, rfca_coord_t pivot ) {
    printf( "Masking: " );
    for( int i =0; i < p->size; i++ ) {
        // For each offset, compute its location on the automaton
        rfca_coord_t c = pattern_offset_abs( pivot, p->offsets[i] );

        rfca_node_t value = rfca_valueC( r, c ) | MASKED_VALUE;
        //if( !(value & MASKED_VALUE) ) {
            rfca_setValueC( r, c, value );
       // }
        printf( "(%d,%d), ", c.row, c.col );
    }
    printf( "\n" );
}

void
unmaskAll( rfca_t* r ) {
    for( int i =0; i < r->buffer->rowCount; i++ ) {
        for( int j =0; j < rfca_buffer_rowLength( r->buffer, i ); j++ ) {
            rfca_node_t value = rfca_buffer_value( r->buffer, i, j );
            if( value & MASKED_VALUE ) {
                rfca_buffer_setValue( r->buffer, i, j, value & ~MASKED_VALUE );
            }
        }
    }
}

int
computeUsage( encoded_rfca_t* v, pattern_t* p ) {
    rfca_t* r = v->rfca;
    pattern_bounds_t pb = pattern_computeBounds( p );
    int usage =0;
    
    for( int i =-pb.rowMin; i < r->buffer->rowCount - pb.rowMax; i++ ) {
        for( int j =-pb.colMin; j < rfca_rowLength( r, i ) - pb.colMax; j++ ) {
            rfca_coord_t pivot = { i,j };
            if( pattern_isMatch( p, r, pivot, 0 ) ) {
                maskRegion( r, p, pivot );
                usage++;
            }
        }
    }
    unmaskAll( r );
    return usage;
}

encoded_rfca_t*
encoded_create_from( rfca_t* r ) {
    // We're creating an encoded version of r using a standard code table
    encoded_rfca_t* v = (encoded_rfca_t*)malloc( sizeof( encoded_rfca_t ) );
    v->rfca =r;

    // Create a standard code table first
    v->codeTable = standardCodeTable( r->opts.base );

    // A lookup table is easier and faster than a list
    pattern_t* pattern_lookup[ r->opts.base+1 ];
    int i =0;
    struct list_head* pos;
    list_for_each( pos, &(v->codeTable->list) ) {
        pattern_list_t* tmp = list_entry( pos, pattern_list_t, list );
        pattern_lookup[i++] = tmp->pattern;
    }

    // The encoded data is represented in a linked list
    v->encoded = (region_list_t*)malloc( sizeof( region_list_t ) );
    INIT_LIST_HEAD( &(v->encoded->list) );

    // However, we need structural information that is not encoded in this list
    // We use a rfca_buffer object to store pointers to the respective regions
    v->index = rfca_buffer_create( r->buffer->width, r->opts.mode );

    // Now we encode each node in the automaton using the standard code table
    for( int i =0; i < r->buffer->rowCount; i++ ) {
        rfca_row_t* row = &r->buffer->rows[i];
        for( int j =0; j < row->size; j++ ) {
            region_list_t* tmp;
            tmp = (region_list_t*)malloc( sizeof( region_list_t ) );
            INIT_LIST_HEAD( &(tmp->list) );

            // Create a region for every singleton on every node
            int value = rfca_value( r, i, j );
            rfca_coord_t pivot = { i,j };
            tmp->region = region_create( pattern_lookup[value], pivot );

            // Add to the encoded dataset
            list_add_tail( &(tmp->list), &(v->encoded->list) );

            // Also add to the indexed dataset
            // Note that we abuse the rfca_buffer's 'value' field by putting a pointer in it
            rfca_buffer_setValueC( v->index, pivot, (uint64_t)tmp->region );

            // Increment the pattern's usage so we can compute its code length later
            pattern_lookup[value]->usage++;
        }
    }

    // Compute the initial encoding sizes for the data and the code table
    v->stdBitsPerSingleton = -log2( 1.0 / (double)r->opts.base );
    pattern_list_updateCodeLength( v->codeTable, v->index->nodeCount );
    v->ctBits = computeCodeTableBits( v );
    v->encodedBits = computeEncodedBits( v );
    return v;
}

void
encoded_free( encoded_rfca_t* v ) {
    region_list_free( v->encoded, true );
    pattern_list_free( v->codeTable, true );
    rfca_buffer_free( v->index );
    free( v );
}

int
encoded_step( encoded_rfca_t* v ) {
    // Test piece
    
    pattern_t* p_union =NULL;
    region_t* r1 = (region_t*)rfca_buffer_value( v->index, 0, 0 );
    region_t* r2 = (region_t*)rfca_buffer_value( v->index, 0, 1 );

    pattern_offset_t p2_offset = pattern_offset( r1->pivot, r2->pivot );
    p_union =pattern_createUnion( r1->pattern, r2->pattern, p2_offset );

    int usage =computeUsage( v, p_union );
    printf( "p_union usage: %d\n", usage );
    
    return 0;
}


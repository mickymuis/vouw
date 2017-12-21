/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include "encoded_rfca.h"
#include <stdlib.h>
#include <math.h>

pattern_list_t*
standardCodeTable( int base ) {
    pattern_list_t* ct;
    ct = (pattern_list_t*)malloc( sizeof( pattern_list_t ) );
    INIT_LIST_HEAD( &(ct->list) );

    for( int i =0; i < base; i++ ) {
        pattern_list_t* tmp;
        tmp = (pattern_list_t*)malloc( sizeof( pattern_list_t ) );
        INIT_LIST_HEAD( &(tmp->list) );

        tmp->pattern = pattern_create_single( i );
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
    pattern_list_computeCodeLength( v->codeTable, v->index->nodeCount );
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

}


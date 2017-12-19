/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include "encoded_rfca.h"
#include <stdlib.h>

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

    // Now we encode each node in the automaton using the standard code table
    for( int i =0; i < r->rowCount; i++ ) {
        rfca_row_t* row = &r->rows[i];
        for( int j =0; j < row->size; j++ ) {
            region_list_t* tmp;
            tmp = (region_list_t*)malloc( sizeof( region_list_t ) );
            INIT_LIST_HEAD( &(tmp->list) );

            int value = rfca_value( r, i, j );
            rfca_coord_t pivot = { i,j };

            tmp->region = region_create( pattern_lookup[value], pivot );
            list_add_tail( &(tmp->list), &(v->encoded->list) );
        }
    }
    return v;
}

void
encoded_free( encoded_rfca_t* v ) {
    region_list_free( v->encoded, true );
    pattern_list_free( v->codeTable, true );
    free( v );
}


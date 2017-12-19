/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include "module_encode.h"
#include "encoded_rfca.h"
#include "list.h"
#include <stdio.h>

void
encoded_print( encoded_rfca_t* v  ) {
    // Label all the patterns so we can print them
    pattern_list_setLabels( v->codeTable );

    // Let's abuse the RFCA datastructure to hold our rendered structure
    rfca_opts_t opts = v->rfca->opts;
    opts.right =true;
    rfca_t* print = rfca_create( opts );

    // Expand each encoded region on the automaton's structure
    struct list_head* pos;
    list_for_each( pos, &(v->encoded->list) ) {
        region_list_t* entry = list_entry( pos, region_list_t, list );
        region_t* region = entry->region;
        pattern_t* pattern = region->pattern;

        for( int i =0; i < pattern->size; i++ ) {
            rfca_coord_t coord = pattern_offset_abs( region->pivot, pattern->offsets[i] );
            rfca_coord_setValue( print, coord, (uint8_t)pattern->label );
        }
    }

    // Pretty print the result
    for( int i =0; i < print->rowCount; i++ ) {
        for( int k =0; k < i * (print->opts.mode-1); k++ )
            fprintf( stdout, " " );
        for( int j =0; j < rfca_rowLength( print, i ); j++ ) {
            fprintf( stdout, "%c", (char)rfca_value( print, i, j ) );
            fprintf( stdout, " " );
        }
        fprintf( stdout, "\n" );
    }

    rfca_free( print );
}

int module_encode( rfca_opts_t opts, int argc, char** argv ) {
    rfca_t* r = rfca_create( opts );
    rfca_generate( r );

    encoded_rfca_t* v = encoded_create_from( r );
    encoded_print( v );
    
    rfca_free( r );
    encoded_free( v );
}


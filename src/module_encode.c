/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include "module_encode.h"
#include "module_print.h"
#include "encoded_rfca.h"
#include "list.h"
#include <stdio.h>

void
encoded_print( encoded_rfca_t* v  ) {
    // Label all the patterns so we can print them
    pattern_list_setLabels( v->codeTable );


    // Let's abuse the RFCA buffer to hold our rendered structure
    rfca_buffer_t* print = rfca_buffer_create( v->rfca->buffer->width, v->rfca->opts.mode  );

    // Expand each encoded region on the automaton's structure
    struct list_head* pos;
    list_for_each( pos, &(v->encoded->list) ) {
        region_t* region = list_entry( pos, region_t, list );
        pattern_t* pattern = region->pattern;

        pattern_setBufferValues( pattern, region->pivot, print, '.' );
        rfca_buffer_setValue( print, region->pivot, pattern->label );
    }

    // Pretty print the result
    /*for( int i =0; i < v->index->rowCount; i++ ) {
        for( int k =0; k < i * (v->index->mode-1); k++ )
            fprintf( stdout, " " );
        for( int j =0; j < rfca_buffer_rowLength( v->index, i ); j++ ) {
            rfca_coord_t c = {i, j};
            region_t* region = (region_t*)rfca_buffer_value( v->index, c );
            fprintf( stdout, "%c", region->pattern->label );
            fprintf( stdout, " " );
        }
        fprintf( stdout, "\n" );
    }*/
    for( int i =0; i < print->rowCount; i++ ) {
        for( int k =0; k < i * (print->mode-1); k++ )
            fprintf( stdout, " " );
        for( int j =0; j < rfca_buffer_rowLength( print, i ); j++ ) {
            rfca_coord_t c = {i, j};
            char label = (char)rfca_buffer_value( print, c );
            fprintf( stdout, "%c", label );
            fprintf( stdout, " " );
        }
        fprintf( stdout, "\n" );
    }
    rfca_buffer_free( print );

    // Print some stats
    fprintf( stdout, "Encoded size: %f bits, code table size: %f bits.\n", 
            v->encodedBits, v->ctBits );

}

void
encoded_printCodeTable( encoded_rfca_t* v ) {
    printf( "Pattern\tSize\tUsage\tCode length\n" );
    printf( "---\t---\t---\t---\n" );
    struct list_head* pos;
    list_for_each( pos, &(v->codeTable->list) ) {
        pattern_t* p= list_entry( pos, pattern_t, list );
        printf( "%c\t%d\t%d\t%f\n",
                p->label, p->size, p->usage, p->codeLength );
    } 
}

int module_encode( rfca_opts_t opts, int argc, char** argv ) {
    rfca_t* r = rfca_create( opts );
    rfca_generate( r );

    encoded_rfca_t* v = encoded_createFrom( r );
    double uncompressed = v->ctBits + v->encodedBits;
    encoded_print( v );
   
   // encoded_test( v );
   // encoded_step( v );
    /*encoded_step( v );
    encoded_step( v );*/
    encoded_encode( v );
    encoded_print( v );
    double compressed = v->ctBits + v->encodedBits;
    printf( "Compression ratio: %f%%\n", compressed / uncompressed * 100.0 );
    encoded_printCodeTable( v );

    rfca_t* r_prime = encoded_decode( v );
    rfca_print( r_prime, true );
    printf( "Correct output? %s\n", rfca_buffer_isEqual( r->buffer, r_prime->buffer ) ? "yes" : "no" );

    rfca_free( r );
    encoded_free( v );
}


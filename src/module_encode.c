/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include "module_encode.h"
#include "module_print.h"
#include "vouw.h"
#include "list.h"
#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "cli.h"

void
vouw_print( vouw_t* v  ) {
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
vouw_printCodeTable( vouw_t* v ) {
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
    rfca_t* r1 = rfca_create( opts );
    rfca_generate( r1 );

    rfca_t* r2 = r1;
    rfca_opts_t opts2 = opts;

    if( argc > 0 && strcmp( argv[0], "using" ) == 0 ) {

        argv++; argc--;
        if( !cli_parseOpts( &opts2, &argv, &argc ) ) {
            rfca_free( r1 );
            return -1;
        }
        r2 = rfca_create( opts2 );
        rfca_generate( r2 );
    }
    
    fprintf( stderr, "RFCA:  %d.%d.%"PRIu64" (%d fold)\n",
            opts2.mode, opts2.base, opts2.rule, opts2.folds );

    vouw_t* v = vouw_createFrom( r2 );
    double uncompressed = v->ctBits + v->encodedBits;
    vouw_encode( v );
    vouw_print( v );
    double compressed = v->ctBits + v->encodedBits;
    printf( "Compression ratio: %f%%\n", compressed / uncompressed * 100.0 );
    vouw_printCodeTable( v );

    rfca_t* r_prime = vouw_decode( v );
    rfca_print( r_prime, true );
    printf( "Correct output? %s\n", rfca_buffer_isEqual( r2->buffer, r_prime->buffer ) ? "yes" : "no" );

    rfca_free( r_prime );

    if( r1 != r2 ) {
        fprintf( stderr, "Now encoding RFCA:  %d.%d.%"PRIu64" (%d) using RFCA: %d.%d.%"PRIu64" (%d)\n",
            opts.mode, opts.base, opts.rule, opts.folds,
            opts2.mode, opts2.base, opts2.rule, opts2.folds);

        vouw_t* v2 = vouw_createEncodedUsing( r1, v->codeTable );
        vouw_print( v2 );
        double compressed = v2->ctBits + v2->encodedBits;
        printf( "Compression ratio: %f%%\n", compressed / uncompressed * 100.0 );
        
        rfca_t* r_prime = vouw_decode( v2 );
        rfca_print( r_prime, true );
        printf( "Correct output? %s\n", rfca_buffer_isEqual( r1->buffer, r_prime->buffer ) ? "yes" : "no" );

        rfca_free( r2 );
        rfca_free( r_prime );
        vouw_free( v2 );
    }
    
    rfca_free( r1 );
    vouw_free( v );
}


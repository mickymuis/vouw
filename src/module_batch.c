/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include "module_batch.h"
#include "vouw.h"
#include "list.h"
#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "cli.h"


int module_encodeAll( rfca_opts_t opts, int argc, char** argv ) {
    uint64_t rulespace = rfca_maxRules( opts.base, opts.mode );
    
    rfca_opts_t opts2 = opts;
    vouw_t* using = NULL;
    double using_baseline =0.0;

    if( argc > 0 && strcmp( argv[0], "using" ) == 0 ) {

        rfca_t* r2 = NULL;
        argv++; argc--;
        if( !cli_parseOpts( &opts2, &argv, &argc ) ) {
            return -1;
        }
        r2 = rfca_create( opts2 );
        rfca_generate( r2 );
        using = vouw_createFrom( r2 );
        vouw_encode( using );

        // Create a baseline by compressing this rule by its own code table
        vouw_t* v2 =vouw_createEncodedUsing( r2, using->codeTable );
        using_baseline = v2->ctBits + v2->encodedBits;
        vouw_free( v2 );

        rfca_free( r2 );
        fprintf( stderr, "Now encoding RFCA class: %d.%d for %"PRIu64" rules, using RFCA: %d.%d.%"PRIu64"\n",
            opts.mode, opts.base, opts.rule,
            opts2.mode, opts2.base, opts2.rule);
    } 
    else
        fprintf( stderr, "Now encoding RFCA class: %d.%d for %"PRIu64" rules\n",
            opts.mode, opts.base, rulespace );

    for( uint64_t i=0; i < rulespace; i++ ) {
        opts.rule = i;
        rfca_t* r =rfca_create( opts );
        rfca_generate( r );

        fprintf( stderr, "Encoding %"PRIu64" (%1.f%%)...", i, (double)i/(double)rulespace * 100.0 );

        vouw_t* v =vouw_createFrom( r );
        double uncompressed = v->ctBits + v->encodedBits;
        vouw_encode( v );
        double compressed = v->ctBits + v->encodedBits;
        double compressed_using =0.0;
        
        if( using ) {
            vouw_free( v );
            v = vouw_createEncodedUsing( r, using->codeTable );
            compressed_using = v->ctBits + v->encodedBits;
        }
        
        fprintf( stderr, "done.\n" );
        printf( "%"PRIu64" %f%%", i, compressed_using / uncompressed * 100.0 );
        if( using )
            printf( " %f%%\n", (compressed_using - compressed) / compressed * 100.0 );
        else
            printf( "\n" );
        

        vouw_free( v );
        rfca_free( r );
    }
}

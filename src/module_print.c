/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include "module_print.h"
#include "ttable.h"
#include <inttypes.h>
#include <stdio.h>

void
rfca_print( rfca_t* r, bool pretty ) {
    for( int i =0; i < r->buffer->rowCount; i++ ) {
        if( pretty ) {
            for( int k =0; k < i * (r->opts.mode-1); k++ )
                fprintf( stdout, " " );
        }
        for( int j =0; j < rfca_rowLength( r, i ); j++ ) {
            fprintf( stdout, "%d", (int)rfca_value( r, i, j ) );
            if( pretty ) fprintf( stdout, " " );
        }
        fprintf( stdout, "\n" );
    }
}

int 
module_print( rfca_opts_t opts, int argc, char** argv ) {

    rfca_t* r = rfca_create( opts );
    rfca_generate( r );
    rfca_print( r, true );
    rfca_free( r );
    return 0;
}

void
ttable_print( ttable_t* tt ) {
    for( int i =0; i < tt->size; i++ ) {
        fprintf( stdout, "[ " );
        for( int j =0; j < tt->inSize; j++ )
            fprintf( stdout, "%d ", tt->entry[i].in[j] );
        fprintf( stdout, "] -> " );
        for( int j =0; j < tt->outSize; j++ )
            fprintf( stdout, "%d ", tt->entry[i].out[j] );
        fprintf( stdout, "\n" );
    }
}

int module_printTTable( rfca_opts_t opts, int argc, char** argv ) {
    
    fprintf( stdout, "Transition table for #%"PRIu64", base %d, mode %d:\n", opts.rule, opts.base, opts.mode );
    ttable_t* tt= ttable_create( opts.base, opts.mode, opts.rule );
    ttable_print( tt );
    ttable_free( tt );
    return 0;
}

void
ttable_printLevel2( ttable_t* tt ) {
    for( int i =0; i < tt->size; i++ ) {
        fprintf( stdout, "[ " );
        for( int j =0; j < tt->mode; j++ )
            fprintf( stdout, "%d ", tt->entry[i].in[j] );
        fprintf( stdout, "] -> %d\n", tt->entry[i].out[0] );
    }
}

int module_printTTable2( rfca_opts_t opts, int argc, char** argv ) {
    
    fprintf( stdout, "Level 2 transition table for #%"PRIu64", base %d, mode %d:\n", opts.rule, opts.base, opts.mode );
    ttable_t* tt2= ttable_createLevel2( opts.base, opts.mode, opts.rule );
    ttable_print( tt2 );
    ttable_free( tt2 );
    return 0;

}

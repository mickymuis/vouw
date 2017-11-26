/*
 * VOUW - Generation, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include "module_print.h"
#include <stdio.h>

void
rfca_print( rfca_t* r ) {
    for( int i =0; i < r->rowCount; i++ ) {
        for( int j =0; j < rfca_rowLength( r, i ); j++ ) {
            fprintf( stdout, " %d", rfca_value( r, i, j ) );
        }
        fprintf( stdout, "\n" );
    }
}

int 
print_main( rfca_opts_t opts, int argc, char** argv ) {

    printf( "PRINT module:\n" );

    rfca_t* r = rfca_create( opts );
    rfca_generate( r );
    rfca_print( r );
    
}


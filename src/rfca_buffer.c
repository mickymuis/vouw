/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include "rfca_buffer.h"
#include <stdlib.h>
#include <string.h>

rfca_buffer_t*
rfca_buffer_create( int width, int mode ) {
    rfca_buffer_t* b = (rfca_buffer_t*)malloc( sizeof( rfca_buffer_t ) );
    
    // We will preallocate everything, growing/shrinking is NOT support for performance reasons
    // Calculate the final number of rows and columns
    b->width = width;
    b->mode = mode;
    b->rowCount = 1;
    int i =width;
    b->nodeCount = i;
    while( i >= mode ) {
        i -= mode-1;
        b->rowCount++;
        b->nodeCount += i;
    }

    // Allocate and zero all rows
    b->rows = malloc( sizeof( rfca_row_t ) * b->rowCount );
    int rowLength = width;
    for( i =0; i < b->rowCount; i++ ) {
        b->rows[i].size = rowLength;
        b->rows[i].cols = malloc( sizeof( rfca_node_t ) * rowLength );
        memset( b->rows[i].cols, 0, rowLength );
        rowLength -= mode-1;
    }

    return b;
}

void
rfca_buffer_free( rfca_buffer_t* b ) {
    for( int i =0; i < b->rowCount; i++ ) {
        free( b->rows[i].cols );
    }
    free( b->rows );
    free( b );
}


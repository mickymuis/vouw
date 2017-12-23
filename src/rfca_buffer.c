/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include "rfca_buffer.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>

rfca_buffer_t*
rfca_buffer_create( int width, int mode ) {
    rfca_buffer_t* b = (rfca_buffer_t*)malloc( sizeof( rfca_buffer_t ) );
    
    // We will preallocate everything, growing/shrinking is NOT supported for performance reasons
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

/*
 * Returns the value of the node at coordinate c
 */
rfca_node_t 
rfca_buffer_value( const rfca_buffer_t* b, rfca_coord_t c ) {
    assert( c.col >= 0 && c.col < b->rows[c.row].size );
    assert( c.row >= 0 && c.row < b->rowCount );
    return b->rows[c.row].cols[c.col];
}

/*
 * Changes the value on coordinate c to value
 */
void
rfca_buffer_setValue( rfca_buffer_t* b, rfca_coord_t c, rfca_node_t value ) {
    assert( c.col >= 0 && c.col < b->rows[c.row].size );
    assert( c.row >= 0 && c.row < b->rowCount );
    b->rows[c.row].cols[c.col] = value;
}

/*
 * Returns the number of columns in row `row'
 */
int 
rfca_buffer_rowLength( const rfca_buffer_t* b, int row ) {
    if( row >= b->rowCount )
        return 0;
    return b->rows[row].size;
}

/*
 * Return true if the coordinate c is within bounds of b
 */
bool
rfca_buffer_checkBounds( const rfca_buffer_t* b, rfca_coord_t c ) {
    return ( c.row >= 0 && c.row < b->rowCount 
             && c.col >= 0 && c.col < b->rows[c.row].size );
}

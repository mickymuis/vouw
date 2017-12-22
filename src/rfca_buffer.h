/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#ifndef RFCA_BUFFER_H
#define RFCA_BUFFER_H

#include <stdint.h>
#include <stdbool.h>

typedef uint64_t rfca_node_t;

typedef struct {
    int row;
    int col;
} rfca_coord_t;

typedef struct {
    rfca_node_t* cols;
    int size;
} rfca_row_t;

typedef struct {
    rfca_row_t* rows;
    int mode;
    int rowCount;
    int width;
    int nodeCount;
} rfca_buffer_t;

rfca_buffer_t*
rfca_buffer_create( int width, int mode );

void
rfca_buffer_free( rfca_buffer_t* b );

rfca_node_t 
rfca_buffer_value( const rfca_buffer_t* b, int row, int col );

rfca_node_t 
rfca_buffer_valueC( const rfca_buffer_t* b, rfca_coord_t c );

void
rfca_buffer_setValue( rfca_buffer_t* b, int row, int col, rfca_node_t value );

void
rfca_buffer_setValueC( rfca_buffer_t* b, rfca_coord_t c, rfca_node_t value );

int 
rfca_buffer_rowLength( const rfca_buffer_t* b, int row );

bool
rfca_buffer_checkBounds( const rfca_buffer_t* r, int row, int col );

bool
rfca_buffer_checkBoundsC( const rfca_buffer_t* b, rfca_coord_t c );

#endif

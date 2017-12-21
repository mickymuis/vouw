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
rfca_buffer_free( rfca_buffer_t* r );

#endif

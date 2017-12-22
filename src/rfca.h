/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#ifndef RFCA_H
#define RFCA_H

#include "rfca_buffer.h"
#include <stdint.h>
#include <stdbool.h>

typedef struct {
    int base;
    int mode;
    uint64_t rule;
    rfca_node_t *input;
    int inputSize;
    int folds;
    bool right; // right-folding automaton
} rfca_opts_t;

typedef struct {
    rfca_opts_t opts;
    int folds;
    rfca_coord_t cur;
    rfca_node_t* ttable;
    rfca_buffer_t* buffer;
} rfca_t;

rfca_t*
rfca_create( rfca_opts_t opts );

void
rfca_free( rfca_t* r );

uint64_t
rfca_maxRules( int base, int mode );

void
rfca_generate( rfca_t* r );

rfca_node_t 
rfca_value( const rfca_t* r, int row, int col );

rfca_node_t 
rfca_valueC( const rfca_t* r, rfca_coord_t c );

void
rfca_setValue( rfca_t* r, int row, int col, rfca_node_t value );

void
rfca_setValueC( rfca_t* r, rfca_coord_t c, rfca_node_t value );

int 
rfca_rowLength( const rfca_t* r, int row );

bool
rfca_checkBounds( const rfca_t* r, int row, int col );

bool
rfca_checkBoundsC( const rfca_t* r, rfca_coord_t c );

uint64_t
pow64( uint64_t a, uint64_t b ); // TODO move to utility

#endif


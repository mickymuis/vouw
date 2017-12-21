/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#ifndef TTABLE_H
#define TTABLE_H

#include <stdint.h>
#include <stdbool.h>
#include "rfca_buffer.h"

void
varbase_incr( rfca_node_t A[], int base, int length );

// Simple array-based transition tables

rfca_node_t*
tt_make( int base, int mode, uint64_t rule );
int
tt_index( int base, int mode, rfca_node_t* A );

// Elaborate transition tables, level 1

typedef struct {
    rfca_node_t* in;
    rfca_node_t* out;
} ttable_entry_t;

typedef struct {
    ttable_entry_t* entry;
    int size;
    int mode;
    int base;
    uint64_t rule;
    int inSize;
    int outSize;
} ttable_t;

ttable_t*
ttable_create( int base, int mode, uint64_t rule );

void
ttable_free( ttable_t* tt );

// Level 2 transition tables

ttable_t*
ttable_createLevel2( int base, int mode, uint64_t rule );

#endif

/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#ifndef PATTERN_H
#define PATTERN_H

#include "list.h"
#include "rfca.h"
#include <stdbool.h>

typedef struct {
    int row;
    int col;
    int value;
} pattern_offset_t;

rfca_coord_t
pattern_offset_abs( rfca_coord_t c, pattern_offset_t offs );

typedef struct {
    pattern_offset_t* offsets;
//    size_t blockSize;
    unsigned int size;
    char label; // for debug printing
} pattern_t;

typedef struct {
    pattern_t* pattern;
    struct list_head list;
} pattern_list_t;

pattern_t*
pattern_create_single( int value );

void
pattern_free( pattern_t* p );

pattern_list_t*
pattern_list_createHead();

void
pattern_list_free( pattern_list_t* list, bool free_patterns );

void 
pattern_list_setLabels( pattern_list_t* list );


#endif

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
    rfca_node_t value;
} pattern_offset_t;

rfca_coord_t
pattern_offset_abs( rfca_coord_t c, pattern_offset_t offs );

pattern_offset_t
pattern_offset( rfca_coord_t pivot1, rfca_coord_t pivot2 );

typedef struct {
    struct list_head list;
    pattern_offset_t* offsets;
    unsigned int usage;
    unsigned int size;
    double codeLength;
    char label; // for debug printing
} pattern_t;

typedef struct {
    int rowMin;
    int rowMax;
    int colMin;
    int colMax;
} pattern_bounds_t;

pattern_t*
pattern_createSingle( int value );

pattern_t*
pattern_createUnion( const pattern_t* p1, const pattern_t* p2, pattern_offset_t p2_offset );

void
pattern_free( pattern_t* p );

void
pattern_updateCodeLength( pattern_t* p, unsigned int totalNodeCount );

double 
pattern_computeCodeLength( const pattern_t* p, unsigned int totalNodeCount );

bool
pattern_isMatch( const pattern_t* p, const rfca_t* r, rfca_coord_t pivot, int variant );

pattern_bounds_t
pattern_computeBounds( const pattern_t* p );

void
pattern_setBufferValues( const pattern_t* p, rfca_coord_t pivot, rfca_buffer_t* b, rfca_node_t value );

void
pattern_list_free( pattern_t* list );

void 
pattern_list_setLabels( pattern_t* list );

double
pattern_list_updateCodeLength( pattern_t* list, unsigned int totalNodeCount );

void
pattern_list_sortByUsageDesc( pattern_t* head );

#endif

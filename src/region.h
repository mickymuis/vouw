/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#ifndef REGION_H
#define REGION_H

#include "rfca.h"
#include "pattern.h"
#include "list.h"

typedef struct {
    struct list_head list;
    pattern_t* pattern;
    rfca_coord_t pivot;
    int variant;
} region_t;

region_t*
region_create( pattern_t* pattern, rfca_coord_t pivot );

void
region_apply( const region_t* region, rfca_t* r );

void
region_free( region_t* r );

void
region_list_free( region_t* );

#endif


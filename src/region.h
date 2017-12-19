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
    pattern_t* pattern;
    rfca_coord_t pivot;
    int variant;
} region_t;

typedef struct {
    region_t* region;
    struct list_head list;
} region_list_t;

region_t*
region_create( pattern_t* pattern, rfca_coord_t pivot );

void
region_free( region_t* r );

region_t*
region_merge( region_t* dst, region_t* src );

region_list_t*
region_list_create();

void
region_list_free( region_list_t*, bool free_regions );

#endif


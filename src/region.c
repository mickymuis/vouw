/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include "region.h"
#include <stdlib.h>

region_t*
region_create( pattern_t* pattern, rfca_coord_t pivot ) {
    region_t* r = (region_t*)malloc( sizeof( region_t ) );
    r->pivot =pivot;
    r->pattern = pattern;
    return r;
}

void
region_free( region_t* r ) {
    free( r );
}

region_t*
region_merge( region_t* dst, region_t* src ) {

}

region_list_t*
region_list_create() {
    region_list_t* list = (region_list_t*)malloc( sizeof( region_list_t ) );
    INIT_LIST_HEAD( &(list->list ) );
    return list;
}

void
region_list_free( region_list_t* list, bool free_regions ) {
    struct list_head* tmp,* pos;
    list_for_each_safe( pos, tmp, &(list->list) ) {
        region_list_t* entry = list_entry( pos, region_list_t, list );
        list_del( pos );
        if( free_regions );
            free( entry->region );
        free( entry );
    }

}

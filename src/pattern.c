/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include "pattern.h"
#include <stdlib.h>
#include <math.h>
//#define BLOCKSIZE 16

rfca_coord_t
pattern_offset_abs( rfca_coord_t c, pattern_offset_t offs ) {
    rfca_coord_t c2 = { c.row + offs.row, c.col + offs.col };
    return c2;
}

pattern_t*
pattern_create_single( int value ) {
    pattern_t* p = (pattern_t*)malloc( sizeof( pattern_t ) );
    
    p->size =1;
    p->usage =0;
//    p->blockSize = BLOCKSIZE;
    p->offsets = (pattern_offset_t*)malloc( sizeof( pattern_offset_t* ) );
    p->offsets->row =0;
    p->offsets->col =0;
    p->offsets->value = value;

    return p;
}

void
pattern_free( pattern_t* p ) {
    free( p->offsets );
    free( p );
}

void
pattern_updateCodeLength( pattern_t* p, unsigned int totalNodeCount ) {
    double probability = (double)p->usage / (double)totalNodeCount;
    p->codeLength = -log2( probability );
}

pattern_list_t*
pattern_list_createHead() {
    pattern_list_t* list;
    list = (pattern_list_t*)malloc( sizeof( pattern_list_t ) );
    INIT_LIST_HEAD( &(list->list) );
    return list;
}

void
pattern_list_free( pattern_list_t* list, bool free_patterns ) {
    struct list_head* tmp,* pos;
    list_for_each_safe( pos, tmp, &(list->list) ) {
        pattern_list_t* entry = list_entry( pos, pattern_list_t, list );
        list_del( pos );
        if( free_patterns );
            free( entry->pattern );
        free( entry );
    }
}

void 
pattern_list_setLabels( pattern_list_t* list ) {
    struct list_head* pos;
    char c = 'A';
    list_for_each( pos, &(list->list) ) {
        pattern_list_t* entry = list_entry( pos, pattern_list_t, list );
        entry->pattern->label =c;
        c++;
    }
}

double
pattern_list_computeCodeLength( pattern_list_t* list, unsigned int totalNodeCount ) {
    struct list_head* pos;
    double totalCodeLength =0.0;
    list_for_each( pos, &(list->list) ) {
        pattern_list_t* entry = list_entry( pos, pattern_list_t, list );
        pattern_updateCodeLength( entry->pattern, totalNodeCount );
        totalCodeLength += entry->pattern->codeLength;
    }
    return totalCodeLength;
}

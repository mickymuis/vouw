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

pattern_offset_t
pattern_offset( rfca_coord_t pivot1, rfca_coord_t pivot2 ) {
    pattern_offset_t o;
    o.value =0;
    o.row = pivot2.row - pivot1.row;
    o.col = pivot2.col - pivot1.col;
    return o;
}

pattern_t*
pattern_createSingle( int value ) {
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

pattern_t*
pattern_createUnion( const pattern_t* p1, const pattern_t* p2, pattern_offset_t p2_offset ) {
    pattern_t* p = (pattern_t*)malloc( sizeof( pattern_t ) );
    
    p->size =p1->size + p2->size;
    p->usage =0;
    p->offsets = (pattern_offset_t*)malloc( p->size * sizeof( pattern_offset_t* ) );
    int k =0;
    for( int i=0; i < p1->size; i++, k++ )
        p->offsets[k] = p1->offsets[i];
    for( int j=0; j < p2->size; j++, k++ ) {
        p->offsets[k] = p2->offsets[j];
        p->offsets[k].row += p2_offset.row;
        p->offsets[k].col += p2_offset.col;
    }

    return p;
}

void
pattern_free( pattern_t* p ) {
    free( p->offsets );
    free( p );
}

void
pattern_updateCodeLength( pattern_t* p, unsigned int totalNodeCount ) {
    p->codeLength = pattern_computeCodeLength( p, totalNodeCount );
}

double 
pattern_computeCodeLength( const pattern_t* p, unsigned int totalNodeCount ) {
    double probability = (double)p->usage / (double)totalNodeCount;
    return -log2( probability );
}

/*
 * Return whether p matches the values in r at the given pivot.
 * Optionally, a non-zero variant can be specified that is added to the values obtained
 * from r before comparison.
 */
bool
pattern_isMatch( const pattern_t* p, const rfca_t* r, rfca_coord_t pivot, int variant ) {
    const int base = r->opts.base;
    // Iterate over all offsets in the pattern
    for( int i =0; i < p->size; i++ ) {
        // For each offset, compute its location on the automaton
        rfca_coord_t c = pattern_offset_abs( pivot, p->offsets[i] );

        // Check bounds
        if( !rfca_checkBoundsC( r, c ) )
            return false;

        // Check match
        if( p->offsets[i].value != ( rfca_valueC( r, c ) + variant )/* % base */)
            return false;
    }
    return true;
}

/*
 * Returns a pattern_bounds object containing the minimum row and col offset
 * and the maximum row and col offset of all offsets in p.
 */
pattern_bounds_t
pattern_computeBounds( const pattern_t* p ) {
    pattern_bounds_t pb;
    if( p->size < 1 ) return pb;
    pb.rowMin = p->offsets[0].row;
    pb.rowMax = p->offsets[0].row;
    pb.colMin = p->offsets[0].col;
    pb.colMax = p->offsets[0].col;

    for( int i =1; i < p->size; i++ ) {
        pattern_offset_t* o =&(p->offsets[0]);
        if( o->row > pb.rowMax )
            pb.rowMax = o->row;
        else if( o->row < pb.rowMin )
            pb.rowMin = o->row;
        if( o->col > pb.colMax )
            pb.colMax = o->col;
        else if( o->col < pb.colMin )
            pb.colMin = o->col;
    }
    return pb;
}

//
// List functions below
//

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
pattern_list_updateCodeLength( pattern_list_t* list, unsigned int totalNodeCount ) {
    struct list_head* pos;
    double totalCodeLength =0.0;
    list_for_each( pos, &(list->list) ) {
        pattern_list_t* entry = list_entry( pos, pattern_list_t, list );
        pattern_updateCodeLength( entry->pattern, totalNodeCount );
        totalCodeLength += entry->pattern->codeLength;
    }
    return totalCodeLength;
}

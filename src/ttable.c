/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include "ttable.h"
#include "rfca.h"
#include <string.h>
#include <stdlib.h>

void
varbase_incr( rfca_node_t A[], int base, int length ) {
    for( int i =length-1; i >= 0; i-- ) {
        if( (++(A[i])) < base ) 
            break;
        else
            A[i] =0;
    }
}

/* 
 * Produce a transition table given base, mode and rule#
 * Returns a pointer to an array of length (base^mode)
 * Each index corresponds to the enumerated rule in a canonical ordering 
 */
rfca_node_t*
tt_make( int base, int mode, uint64_t rule ) {
    int rulesize = pow64( base, mode );
    rfca_node_t* tt = malloc( sizeof( rfca_node_t ) * rulesize );

    memset( tt, 0, rulesize );

    uint64_t decimal = rule;
    int i = rulesize - 1;
    while( i >= 0 && decimal > 0 ) {
        tt[i] = decimal % base;
        decimal = decimal / base;
        i--;
    }
    return tt;
}

/* 
 * Given base and mode, return the rule number that indexes the transition table
 * In order to compute the rule number, mode bytes are read from A
 * The bytes in A..A+mode-1 are assumed to be smaller than base 
 */
int
tt_index( int base, int mode, rfca_node_t* A ) {
    int mult =1;
    int index =0;
    for( int i = mode-1; i >= 0; i-- ) {
        index += A[i] * mult;
        mult *= base;
    }
    return index;
}

ttable_t*
ttable_create( int base, int mode, uint64_t rule ) {
    // Prepare a buffer to store the input pattern as we increment it
    rfca_node_t A[mode];
    memset( A, 0, mode );

    // Allocate the table entries depending on the mode/base
    ttable_t* tt = malloc( sizeof( ttable_t ) );
    tt->size = pow64( base, mode );
    tt->entry = malloc( sizeof( ttable_entry_t ) * tt->size );
    tt->base =base;
    tt->mode =mode;
    tt->rule =rule;
    tt->inSize = mode;
    tt->outSize = 1;

    // Populate the input part
    for( int i =0; i < tt->size; i++ ) {
        tt->entry[i].in = malloc( sizeof( rfca_node_t ) * mode );
        memcpy( tt->entry[i].in, A, mode );
        varbase_incr( A, base, mode );
    }

    // Populate the output part
    uint64_t decimal = rule;
    int i = tt->size - 1;
    while( i >= 0 && decimal >= 0 ) {
        tt->entry[i].out = malloc( sizeof( rfca_node_t ) );
        tt->entry[i].out[0] = decimal % base;
        decimal = decimal / base;
        i--;
    }
    return tt;

}

void
ttable_free( ttable_t* tt ) {
    for( int i =0; i < tt->size; i++ ) {
        free( tt->entry[i].in );
        free( tt->entry[i].out );
    }
    free( tt->entry );
    free( tt );
}

ttable_t*
ttable_createLevel2( int base, int mode, uint64_t rule ) {
    // Prepare a buffer to store the double input pattern as we increment it
    rfca_node_t A[mode*mode];
    memset( A, 0, mode*mode );

    // We use a simple array-based tt as input
    rfca_node_t* tt = tt_make( base, mode, rule );

    // Allocate the table entries depending on the mode/base
    ttable_t* tt2 = malloc( sizeof( ttable_t ) );
    int rulesize = pow64( base, mode );
    tt2->size = rulesize * rulesize;
    tt2->entry = malloc( sizeof( ttable_entry_t ) * tt2->size );
    tt2->base =base;
    tt2->mode =mode;
    tt2->rule =rule;
    tt2->inSize = mode*mode;
    tt2->outSize = mode;

    for( int i =0; i < tt2->size; i++ ) {
        // Allocate the input pattern and initialize with A
        rfca_node_t* in = malloc( sizeof( rfca_node_t ) * mode*mode );
        memcpy( in, A, mode*mode );

        // Allocate the output pattern and use the level 1 TT to fill it
        rfca_node_t* out = malloc( sizeof( rfca_node_t ) * mode );
        for( int j =0; j < mode; j++ ) {
            out[j] = tt[tt_index( base, mode, &A[j*mode] )];
        }

        tt2->entry[i].in = in;
        tt2->entry[i].out = out;

        varbase_incr( A, base, mode*mode ); // Increment the input pattern
    }

    return tt2;
}


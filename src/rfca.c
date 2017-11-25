/*
 * VOUW - Generation, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include "rfca.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int row;
    int col;
} rfca_coord_t;

/* Calculate a^b with 64-bit integers */
uint64_t
pow64( uint64_t a, uint64_t b ) {
    if( !b )
        return 1;
    else if( b % 2 == 1 )
        return a * pow64( a, b - 1 );
    /* ELSE */
    uint64_t p = pow64( a, b/2 );
    return p*p;
}

/* Produce a transition table given base, mode and rule#
 * Returns a pointer to an array of length (base^mode)
 * Each index corresponds to the enumerated rule in a canonical ordering */
uint8_t*
makeTTable( int base, int mode, uint64_t rule ) {
    int rulesize = pow64( base, mode );
    uint8_t* tt = malloc( sizeof( uint8_t ) * rulesize );

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

/* Given base and mode, return the rule number that indexes the transition table
 * In order to compute the rule number, mode bytes are read from A
 * The bytes in A..A+mode-1 are assumed to be smaller than base */
int
ttIndex( int base, int mode, uint8_t* A ) {
    int mult =1;
    int index =0;
    for( int i = mode-1; i >= 0; i-- ) {
        index += A[i] * mult;
        mult *= base;
    }
    return index;
}

rfca_t*
rfca_create( rfca_opts_t opts ) {
    rfca_t* r = malloc( sizeof( rfca_t ) );

    r->opts = opts;

    // We will preallocate everything, growing/shrinking is NOT support for performance concerns
    // Calculate the final number of rows and columns
    r->width = opts.inputSize + opts.folds;
    r->rowCount = 1;
    int i =r->width;
    r->nodeCount =i;
    while( i >= opts.mode ) {
        i -= opts.mode-1;
        r->rowCount++;
        r->nodeCount += i;
    }

    // Allocate and zero all rows
    r->rows = malloc( sizeof( rfca_row_t ) * r->rowCount );
    int rowLength = r->width;
    for( i =0; i < r->rowCount; i++ ) {
        r->rows[i].size = rowLength;
        r->rows[i].cols = malloc( sizeof( uint8_t ) * rowLength );
        memset( r->rows[i].cols, 0, rowLength );
    }

    // Write the input at the beginning
    // If right == true, we simply reverse the input in order to keep all other code simpler
    for( i =0; i < opts.inputSize; i++ ) {
        r->rows[0].cols[i] = 
            opts.right ? opts.input[(opts.inputSize - 1)- i] : opts.input[i];

        // Sanity check, user input should also be checked elsewhere!
        if( r->rows[0].cols[i] >= opts.base ) 
            r->rows[0].cols[i] = opts.base-1;
    }

    r->folds = 0; // FIXME remove
    r->curRow = 0;
    r->curPos = opts.inputSize-1; // Position at the last input node

    // Finally, populate the transition table based on base, mode and rule number
    r->ttable = makeTTable( opts.base, opts.mode, opts.rule );

    return r;
}

void
rfca_free( rfca_t* r ) {
    for( int i =0; i < r->rowCount; i++ ) {
        free( r->rows[i].cols );
    }
    free( r->rows );
    free( r->ttable );
    free( r );
}

uint64_t
rfca_maxRules( int mode, int base ) {
    uint64_t rulesize = pow64( base, mode );
    return pow64( base, rulesize );
}



/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include "rfca.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "ttable.h"

#define STEP_REDUCE 1
#define STEP_FOLD 2
#define STEP_DONE 0

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


/*
 * Construct a new rfca_t object given the parameters in opts.
 * All memory will be pre-allocated, thus shrinking/growing is not possible after this call.
 * Also computed the transition table for the given base/mode/rule, but the nodes are not
 * yet computed (see rfca_generate())
 * Returns a pointer to a heap-allocated instance of the object.
 */
rfca_t*
rfca_create( rfca_opts_t opts ) {
    rfca_t* r = malloc( sizeof( rfca_t ) );

    r->opts = opts;
    r->buffer = rfca_buffer_create( opts.inputSize + opts.folds, opts.mode );

    // Write the input at the beginning
    // If right == true, we simply reverse the input in order to keep all other code simpler
    for( int i =0; i < opts.inputSize; i++ ) {
        r->buffer->rows[0].cols[i] = 
            opts.right ? opts.input[i] : opts.input[(opts.inputSize - 1)- i];

        // Sanity check, user input should also be checked elsewhere!
        if( r->buffer->rows[0].cols[i] >= opts.base ) 
            r->buffer->rows[0].cols[i] = opts.base-1;
    }

    r->folds = 0; // FIXME remove
    r->cur.row = 0;
    r->cur.col = opts.inputSize-1; // Position at the last input node

    // Finally, populate the transition table based on base, mode and rule number
    r->ttable = tt_make( opts.base, opts.mode, opts.rule );

    return r;
}

/* 
 * Release all memory occupied by the rfca_t object and its underlying structures 
 */
void
rfca_free( rfca_t* r ) {
    free( r->ttable );
    rfca_buffer_free( r->buffer );
    free( r );
}

/* 
 * Compute the total number of rules given a mode and base 
 */
uint64_t
rfca_maxRules( int base, int mode ) {
    uint64_t rulesize = pow64( base, mode );
    return pow64( base, rulesize );
}

/* 
 * Given a {row,col}, calculate the position of the next node
 */
rfca_coord_t
next( rfca_t* r, rfca_coord_t c ) {
    // Length of this row if only the original input was used
    int unfoldedRowLength = r->opts.inputSize - c.row * (r->opts.mode-1);
    // Last row that contains nodes reduced from original input
    int lastInputRow = r->opts.inputSize / (r->opts.mode-1) - 1;
    if( r->opts.inputSize % (r->opts.mode-1) != 0 )
        lastInputRow++;

    if( c.col < unfoldedRowLength && c.row <= lastInputRow ) {
        if( c.col == unfoldedRowLength - 1 ) {
            if( c.row == lastInputRow ) {
                // First pivot
                // { col: inputSize, row: 0 }
                c.col = r->opts.inputSize;
                c.row = 0;
            } else {
                // Next row
                // { col: 0, row: row + 1 }
                c.col = 0;
                c.row++;
            }
        }
        else {
            // Next col
            // { col: col + 1, row: row };
            c.col++;
        }
    }
    else if( c.col == 0 ) {
        // Pivot -> fold
        // { col: inputSize + (row - lastInputRow), row: 0 }
        c.col = r->opts.inputSize + (c.row - lastInputRow);
        c.row = 0;

    } 
    else {
        // Next row in when folding
        // { col: col - 1, row: row + 1 }
        c.col--;
        c.row++;
    }
    return c;
}

/* 
 * Advance the automaton by one step 
 */
int 
step( rfca_t* r ) {
    //fprintf( stdout, "{row: %d, col: %d}\n", r->cur.row, r->cur.col );
    // Step 0.: compute the node position that needs to be updated
    rfca_coord_t nextPos = next( r, r->cur );
    if( nextPos.col >= r->buffer->width || nextPos.row >= r->buffer->rowCount )
        return STEP_DONE;

    // Step 1a. check whether a special action is required
    if( nextPos.row == 0 ) {
        // Top row, which means the next step can only be a fold
        // Step 1b. copy ('fold') the apex/singleton row over to the top row
        int i = r->cur.row;
        rfca_node_t foldValue = r->buffer->rows[i].cols[0];
        r->buffer->rows[0].cols[nextPos.col] = foldValue;
        
        r->folds++;
        r->cur = nextPos;
        return STEP_FOLD;
    }
    // Now we calculate the next iteration by reduction
    // Step 1. get values from parent nodes (one row up)
    rfca_node_t* parents = r->buffer->rows[nextPos.row-1].cols + nextPos.col;

    // Step 2. Use the transition table to obtain the value for the current node
    int i = tt_index( r->opts.base, r->opts.mode, parents );
    rfca_node_t value =r->ttable[i];
    r->buffer->rows[nextPos.row].cols[nextPos.col] = value;

    r->cur = nextPos;
    return STEP_REDUCE;
}

/*
 * Compute the values for all allocated nodes 
 */
void
rfca_generate( rfca_t* r ) {
    while( step( r ) != STEP_DONE );
}

/*
 * Returns the transposed (mirrored) coordinates for {col,row}
 * This translates internal coordinates to logical (abstracted) coordinates
 */
rfca_coord_t
transpose( rfca_t* r, rfca_coord_t c ) {
    assert( c.row < r->buffer->rowCount );
    if( r->opts.right )
        return c; // Not mirrored
    // Mirrored 
    c.col = (r->buffer->rows[c.row].size-1) - c.col;
    return c;
}

/*
 * Returns the value of the node at the logical coordinate c
 */
rfca_node_t 
rfca_value( rfca_t* r, int row, int col ) {
    assert( row < r->buffer->rowCount );
    rfca_coord_t c = { row, col };
    return rfca_coord_value( r, c );
}

/*
 * Returns the value of the node at the logical coordinate c
 */
rfca_node_t 
rfca_coord_value( rfca_t* r, rfca_coord_t c ) {
    c = transpose( r, c );
    assert( c.col < r->buffer->rows[c.row].size );
    assert( c.row < r->buffer->rowCount );
    return r->buffer->rows[c.row].cols[c.col];
}

/*
 * Changes the value on logical coordinate {row,col} to value
 * Note that this is not the normal operation of the automaton
 */
void
rfca_setValue( rfca_t* r, int row, int col, rfca_node_t value ) {
    assert( row < r->buffer->rowCount );
    rfca_coord_t c = { row, col };
    rfca_coord_setValue( r, c, value );
}

/*
 * Changes the value on logical coordinate c to value
 * Note that this is not the normal operation of the automaton
 */
void
rfca_coord_setValue( rfca_t* r, rfca_coord_t c, rfca_node_t value ) {
    c = transpose( r, c );
    assert( c.col < r->buffer->rows[c.row].size );
    assert( c.row < r->buffer->rowCount );
    r->buffer->rows[c.row].cols[c.col] = value;
}

/*
 * Returns the number of columns in row `row'
 */
int 
rfca_rowLength( rfca_t* r, int row ) {
    if( row >= r->buffer->rowCount )
        return 0;
    return r->buffer->rows[row].size;
}

/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include "encoded_rfca.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


static pattern_t*
standardCodeTable( int base ) {
    pattern_t* ct;
    ct = (pattern_t*)malloc( sizeof( pattern_t ) );
    INIT_LIST_HEAD( &(ct->list) );

    for( int i =0; i < base; i++ ) {
        pattern_t* tmp;
        tmp = pattern_createSingle( i );

        list_add_tail( &(tmp->list), &(ct->list ) );
    }
    return ct;
}

static double
computeCodeTableBits( encoded_rfca_t* v ) {
    double bits =0.0;
    struct list_head* pos;
    list_for_each( pos, &(v->codeTable->list) ) {
        pattern_t* tmp = list_entry( pos, pattern_t, list );
        bits += tmp->codeLength + v->stdBitsPerSingleton * tmp->size;
    }
    return bits;
}

static double
computeEncodedBits( encoded_rfca_t* v ) {
    double bits =0.0;
    struct list_head* pos;
    list_for_each( pos, &(v->encoded->list) ) {
        region_t* tmp = list_entry( pos, region_t, list );
        bits += tmp->pattern->codeLength + v->stdBitsPerPivot;
    }
    return bits;
}

static void
maskRegion( rfca_t* r, pattern_t* p, rfca_coord_t pivot ) {
//    printf( "Masking: " );
    for( int i =0; i < p->size; i++ ) {
        // For each offset, compute its location on the automaton
        rfca_coord_t c = pattern_offset_abs( pivot, p->offsets[i] );

        rfca_setMasked( r, c, true );
//        printf( "(%d,%d), ", c.row, c.col );
    }
//    printf( "\n" );
}

static int
computeUsage( encoded_rfca_t* v, pattern_t* p ) {
    rfca_t* r = v->rfca;
    pattern_bounds_t pb = pattern_computeBounds( p );
    int usage =0;
    
    for( int i =-pb.rowMin; i < r->buffer->rowCount - pb.rowMax; i++ ) {
        for( int j =-pb.colMin; j < rfca_rowLength( r, i ) - pb.colMax; j++ ) {
            rfca_coord_t pivot = { i,j };
            if( pattern_isMatch( p, r, pivot, 0 ) ) {
                maskRegion( r, p, pivot );
                usage++;
            }
        }
    }
    rfca_unmaskAll( r );
    return usage;
}

/*
 * Calculate the gain in encoding size if patterns p1 and p2 were replaced by their union p.
 * p is assumed to be the union of patterns p1 and p2 and to have estimated usage p_usage.
 * The usage of p is assumed to be less or equal than the usages of p1 and p2.
 * The return value is the difference in encoding side in bits.
 */
static double
computeGain( encoded_rfca_t* v, pattern_t* p1, pattern_t* p2, pattern_t* p, int p_usage ) {

    const int totalNodes = v->rfca->buffer->nodeCount;
    const double oldBits = v->ctBits + v->encodedBits; // MDL's L(M) + L(M|D)
    double newBits = oldBits;

    // Step 1. remove the size of the old patterns completely
    newBits -= (p1->codeLength + v->stdBitsPerPivot) * p1->usage;
    newBits -= (p2->codeLength + v->stdBitsPerPivot) * p2->usage;

    // Step 2. add the bits from the old patterns that still have usage > 0
    int p1_usage = p1->usage - p_usage;
    if( p1_usage > 0 )
        // Pattern p1 will remain in the code table
        newBits += (-log2( (double)p1_usage / (double)totalNodes ) + v->stdBitsPerPivot) * p1_usage;
    else if( p1->size > 1 ) 
        // Pattern p1 will be removed from the code table
        newBits -= (v->stdBitsPerSingleton * p1->size + p1->codeLength);

    int p2_usage = p2->usage - p_usage;
    if( p2_usage > 0 )
        newBits += (-log2( (double)p2_usage / (double)totalNodes ) + v->stdBitsPerPivot) * p2_usage;
    else if( p2->size > 1 ) 
        // Pattern p2 will be removed from the code table
        newBits -= (v->stdBitsPerSingleton * p2->size + p2->codeLength);

    // Step 3. add the length from the union pattern p
    double p_codeLength = -log2( (double)p_usage / (double)totalNodes );  
    // data part
    newBits += (p_codeLength + v->stdBitsPerPivot) * (p_usage);
    // code table part
    newBits += v->stdBitsPerSingleton * p->size + p_codeLength;
    
    return oldBits - newBits;
}

static region_t*
createRegion( encoded_rfca_t* v, pattern_t* p, rfca_coord_t pivot, bool mask ) {
    rfca_t* r = v->rfca;
    region_t* region = (region_t*)malloc( sizeof( region_t ) );
    region->pivot =pivot;
    region->pattern =p;

    for( int i =0; i < p->size; i++ ) {
        // For each offset, compute its location on the automaton
        rfca_coord_t c = pattern_offset_abs( pivot, p->offsets[i] );

        // There is probably an existing region at this location
        region_t* ro = (region_t*)rfca_buffer_value( v->index, c );
        if( ro ) {
            // Substract the usage from the region's pattern
            ro->pattern->usage--;
            // Remove this region from the index
            pattern_setBufferValues( ro->pattern, ro->pivot, v->index, 0 );
            // Remove this region from the list and free its memory
            list_del( &(ro->list) );
            free( ro );
        }

        rfca_buffer_setValue( v->index, c, (uint64_t)region );

        if( mask )
            rfca_setMasked( r, c, true );
    }

    list_add( &(region->list), &(v->encoded->list) );
    return region;
}

/* 
 * Reencode the rfca encoded in v by adding pattern p to the code table
 */
static void
updateEncoding( encoded_rfca_t* v, pattern_t* p ) {
    rfca_t* r = v->rfca;
    pattern_bounds_t pb = pattern_computeBounds( p );
    
    for( int i =-pb.rowMin; i < r->buffer->rowCount - pb.rowMax; i++ ) {
        for( int j =-pb.colMin; j < rfca_rowLength( r, i ) - pb.colMax; j++ ) {
            rfca_coord_t pivot = { i,j };
            if( pattern_isMatch( p, r, pivot, 0 ) ) {
                createRegion( v, p, pivot, true );
                p->usage++;
            }
        }
    }
    rfca_unmaskAll( r );

    list_add_tail( &(p->list), &(v->codeTable->list) );
}

encoded_rfca_t*
encoded_create_from( rfca_t* r ) {
    // We're creating an encoded version of r using a standard code table
    encoded_rfca_t* v = (encoded_rfca_t*)malloc( sizeof( encoded_rfca_t ) );
    v->rfca =r;

    // Create a standard code table first
    v->codeTable = standardCodeTable( r->opts.base );

    // A lookup table is easier and faster than a list
    pattern_t* pattern_lookup[ r->opts.base+1 ];
    int i =0;
    struct list_head* pos;
    list_for_each( pos, &(v->codeTable->list) ) {
        pattern_t* tmp = list_entry( pos, pattern_t, list );
        pattern_lookup[i++] = tmp;
    }

    // The encoded data is represented in a linked list
    v->encoded = (region_t*)malloc( sizeof( region_t ) );
    INIT_LIST_HEAD( &(v->encoded->list) );

    // However, we need structural information that is not encoded in this list
    // We use a rfca_buffer object to store pointers to the respective regions
    v->index = rfca_buffer_create( r->buffer->width, r->opts.mode );

    // Now we encode each node in the automaton using the standard code table
    for( int i =0; i < r->buffer->rowCount; i++ ) {
        rfca_row_t* row = &r->buffer->rows[i];
        for( int j =0; j < row->size; j++ ) {

            // Create a region for every singleton on every node
            rfca_coord_t pivot = { i,j };
            int value = rfca_value( r, pivot );
            region_t* region = region_create( pattern_lookup[value], pivot );

            // Add to the encoded dataset
            list_add_tail( &(region->list), &(v->encoded->list) );

            // Also add to the indexed dataset
            // Note that we abuse the rfca_buffer's 'value' field by putting a pointer in it
            rfca_buffer_setValue( v->index, pivot, (uint64_t)region );

            // Increment the pattern's usage so we can compute its code length later
            pattern_lookup[value]->usage++;
        }
    }

    // Compute the initial encoding sizes for the data and the code table
    v->stdBitsPerSingleton = -log2( 1.0 / (double)r->opts.base );
    v->stdBitsPerPivot = -log2( 1.0 / (double)v->index->nodeCount );
    pattern_list_updateCodeLength( v->codeTable, v->index->nodeCount );
    v->ctBits = computeCodeTableBits( v );
    v->encodedBits = computeEncodedBits( v );
    return v;
}

void
encoded_free( encoded_rfca_t* v ) {
    region_list_free( v->encoded );
    pattern_list_free( v->codeTable );
    rfca_buffer_free( v->index );
    free( v );
}

int
encoded_step( encoded_rfca_t* v ) {

    // Take the two patterns with the largest usage
    pattern_t* p1 = list_entry( v->codeTable->list.next, pattern_t, list );
    pattern_t* p2 = list_entry( p1->list.next, pattern_t, list );

    printf( "usages | p1: %d, p2: %d\n", p1->usage, p2->usage );

    double bestGain =0.0;
    pattern_t* bestP =NULL;
    int bestUsage =0;

    struct list_head *pos1, *pos2;
    list_for_each( pos1, &(v->encoded->list) ) {
        region_t* r1 = list_entry( pos1, region_t, list );
        if( r1->pattern != p1 )
            continue;
        list_for_each( pos2, &(v->encoded->list) ) { 
            region_t* r2 = list_entry( pos2, region_t, list );
            if( r2->pattern != p2 )
                continue;
            
            pattern_offset_t p2_offset = pattern_offset( r1->pivot, r2->pivot );
            pattern_t* p_union =pattern_createUnion( r1->pattern, r2->pattern, p2_offset );
            int usage =computeUsage( v, p_union );
            double gain = computeGain( v, r1->pattern, r2->pattern, p_union, usage );
            if( gain > bestGain ) {
                bestGain =gain;
                if( bestP )
                    pattern_free( bestP );
                bestP = p_union;
                bestUsage =usage;
            }

        }
    }
    if( bestGain == 0.0 ) {
        printf( "encoded_step(): No compression gain.\n" );
        return false;
    }

    printf( "encoded_step(): bestP usage: %d\n", bestUsage );
    printf( "encoded_step(): compression size gain: %f bits\n", bestGain );

    updateEncoding( v, bestP );

    return true;
}

int
encoded_test( encoded_rfca_t* v ) {
    // Test piece
    
    pattern_t* p_union =NULL;
    rfca_coord_t c = {0,0};
    region_t* r1 = (region_t*)rfca_buffer_value( v->index, c );
    c.col = 1;
    region_t* r2 = (region_t*)rfca_buffer_value( v->index, c );

    printf( "Code lengths: p1: %f, p2 %f\n", r1->pattern->codeLength, r2->pattern->codeLength );
    printf( "Total nodes: %d\n", v->index->nodeCount );

    pattern_offset_t p2_offset = pattern_offset( r1->pivot, r2->pivot );
    p_union =pattern_createUnion( r1->pattern, r2->pattern, p2_offset );

    int usage =computeUsage( v, p_union );
    printf( "p_union usage: %d\n", usage );
    double gain = computeGain( v, r1->pattern, r2->pattern, p_union, usage );
    printf( "Compression size gain: %f bits\n", gain );

    updateEncoding( v, p_union );
    pattern_list_updateCodeLength( v->codeTable, v->index->nodeCount );
    v->ctBits = computeCodeTableBits( v );
    v->encodedBits = computeEncodedBits( v );
    
    return 0;
}


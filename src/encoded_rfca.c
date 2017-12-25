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
        bits += tmp->pattern->codeLength + v->stdBitsPerPivot + v->stdBitsPerVariant;
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
computeUsage( encoded_rfca_t* v, pattern_t* p1, int v1, pattern_t* p2, int v2, pattern_offset_t p2_offset ) {
    int usage =0;
    const int base = v->rfca->opts.base;
    struct list_head *pos1, *pos2;

    list_for_each( pos1, &(v->encoded->list) ) {
        region_t* r1 = list_entry( pos1, region_t, list );
        if( r1->pattern != p1 )
            continue;
        // Compute the difference between r1's variant and variant v1
        int vv1 = (v1 - r1->variant) % base;
        vv1 = vv1 < 0 ? vv1+base : vv1;

        list_for_each( pos2, &(v->encoded->list) ) { 
            region_t* r2 = list_entry( pos2, region_t, list );
            if( r2->pattern != p2 )
                continue;
            
            // Compute the difference between r2's variant and variant v2
            int vv2 = (v2 - r2->variant) % base;
            vv2 = vv2 < 0 ? vv2+base : vv2;
            
            /*pattern_offset_t candidate_p2_offset = pattern_offset( r1->pivot, r2->pivot );
            if( candidate_p2_offset.row == p2_offset.row &&
                candidate_p2_offset.col == p2_offset.col )
                usage ++;*/

            if( vv1 == vv2 &&
                r2->pivot.row - r1->pivot.row == p2_offset.row &&
                r2->pivot.col - r1->pivot.col == p2_offset.col )
                usage ++;
        }
    }
    return usage;
}


static int 
computeUsage2( encoded_rfca_t* v, pattern_t* p1, pattern_t* p2, pattern_offset_t p2_offset ) {
    int usage =0;
    struct list_head *pos1, *pos2;
    list_for_each( pos1, &(v->encoded->list) ) {
        region_t* r1 = list_entry( pos1, region_t, list );
        if( r1->pattern != p1 )
            continue;
        list_for_each( pos2, &(v->encoded->list) ) { 
            region_t* r2 = list_entry( pos2, region_t, list );
            if( r2->pattern != p2 || (r1->pivot.row == r2->pivot.row && r1->pivot.col == r2->pivot.col) )
                continue;
            
            /*pattern_offset_t candidate_p2_offset = pattern_offset( r1->pivot, r2->pivot );
            if( candidate_p2_offset.row == p2_offset.row &&
                candidate_p2_offset.col == p2_offset.col )
                usage ++;*/

            if( r2->pivot.row - r1->pivot.row == p2_offset.row &&
                r2->pivot.col - r1->pivot.col == p2_offset.col )
                usage ++;
        }
    }
    return usage;
}

/*
 * Calculate the gain in encoding size if patterns p1 and p2 were replaced by their union.
 * This union is assumed to have estimated usage p_usage,
 * which is assumed to be less or equal than the usages of p1 and p2.
 * The return value is the difference in encoding side in bits.
 */
static double
computeGain( encoded_rfca_t* v, pattern_t* p1, pattern_t* p2, int p_usage ) {

    const int totalNodes = v->rfca->buffer->nodeCount;
    const double oldBits = v->ctBits + v->encodedBits; // MDL's L(M) + L(M|D)
    double newBits = oldBits;

    // Step 1. remove the size of the old patterns completely
    newBits -= (p1->codeLength + v->stdBitsPerPivot + v->stdBitsPerVariant) * p1->usage;
    newBits -= (p2->codeLength + v->stdBitsPerPivot + v->stdBitsPerVariant) * p2->usage;

    // Step 2. add the bits from the old patterns that still have usage > 0
    int p1_usage = p1->usage - p_usage;
    if( p1_usage > 0 )
        // Pattern p1 will remain in the code table
        newBits += (-log2( (double)p1_usage / (double)totalNodes ) + v->stdBitsPerPivot + v->stdBitsPerVariant) * p1_usage;
    else if( p1->size > 1 ) 
        // Pattern p1 will be removed from the code table
        newBits -= (v->stdBitsPerSingleton * p1->size + p1->codeLength);

    int p2_usage = p2->usage - p_usage;
    if( p2_usage > 0 )
        newBits += (-log2( (double)p2_usage / (double)totalNodes ) + v->stdBitsPerPivot + v->stdBitsPerVariant) * p2_usage;
    else if( p2->size > 1 ) 
        // Pattern p2 will be removed from the code table
        newBits -= (v->stdBitsPerSingleton * p2->size + p2->codeLength);

    // Step 3. add the length from the union pattern p
    double p_codeLength = -log2( (double)p_usage / (double)totalNodes );  
    // data part
    newBits += (p_codeLength + v->stdBitsPerPivot + v->stdBitsPerVariant) * (p_usage);
    // code table part
    newBits += v->stdBitsPerSingleton * (p1->size + p2->size) + p_codeLength;
    
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

static pattern_t*
mergeEncodedPatterns( encoded_rfca_t* v, pattern_t* p1, pattern_t* p2, pattern_offset_t p2_offset ) {
    // Create the union pattern of p1 and p2
    pattern_t* p_union = pattern_createUnion( p1, p2, p2_offset );
    list_add( &(p_union->list), &(p2->list) );
    //list_add( &(p_union->list), &(v->codeTable->list) );

    // Iterate over every possible combination of p1 and p2
    // Complexity is approx. (N/2)^2 in the list of regions. I'm not proud of it.
    struct list_head *pos1, *pos2, *tmp1, *tmp2;
    list_for_each_safe( pos1, tmp1, &(v->encoded->list) ) {
        region_t* r1 = list_entry( pos1, region_t, list );
        if( r1->pattern != p1 )
            continue;
        list_for_each_safe( pos2, tmp2, &(v->encoded->list) ) { 
            region_t* r2 = list_entry( pos2, region_t, list );
            if( r2->pattern != p2 )
                continue;
            
            pattern_offset_t candidate_p2_offset = pattern_offset( r1->pivot, r2->pivot );
            if( candidate_p2_offset.row == p2_offset.row &&
                candidate_p2_offset.col == p2_offset.col ) {

                rfca_coord_t pivot = r1->pivot;
                // Fix the list iterators if we're removing r1 and r2
                if( tmp1 == &r2->list )
                    tmp1 = r2->list.next;
                if( tmp2 == &r1->list )
                    tmp2 = r1->list.next;
                // Remove and free both r1 and r2
                r1->pattern->usage--;
                list_del( &(r1->list) );
                region_free( r1 );
                r2->pattern->usage--;
                list_del( &(r2->list) );
                region_free( r2 );
                // Create a new region at this pivot containing p_union
                region_t* region = (region_t*)malloc( sizeof( region_t ) );
                region->pivot =pivot;
                region->pattern =p_union;

                list_add_tail( &(region->list), &(v->encoded->list) );

                p_union->usage ++;
            }
        }
    }

    // If p1 and/or p2 are not used and they are not singleton patterns,
    // we remove them from the code table and free their memory.
    if( p1->usage == 0 && p1->size > 1 ) {
        list_del( &(p1->list ) );
        pattern_free( p1 );
    }
    if( p2->usage == 0 && p2->size > 1 ) {
        list_del( &(p2->list ) );
        pattern_free( p2 );
    }

    return p_union;
}


///
/// 2
///
static pattern_t*
mergeEncodedPatterns2( encoded_rfca_t* v, pattern_t* p1, int v1, pattern_t* p2, int v2, pattern_offset_t p2_offset ) {
    const int base = v->rfca->opts.base;
    // Create the union pattern of p1 and p2
    pattern_t* p_union = pattern_createVariantUnion( p1, v1, p2, v2, p2_offset, base );
    list_add( &(p_union->list), &(p2->list) );
    //list_add( &(p_union->list), &(v->codeTable->list) );

    // Iterate over every possible combination of p1 and p2
    // Complexity is approx. (N/2)^2 in the list of regions. I'm not proud of it.
    struct list_head *pos1, *pos2, *tmp1, *tmp2;
    list_for_each_safe( pos1, tmp1, &(v->encoded->list) ) {
        region_t* r1 = list_entry( pos1, region_t, list );
        if( r1->pattern != p1 )
            continue;
        // Compute the difference between r1's variant and variant v1
        int vv1 = (v1 - r1->variant) % base;
        vv1 = vv1 < 0 ? vv1+base : vv1;
        
        list_for_each_safe( pos2, tmp2, &(v->encoded->list) ) { 
            region_t* r2 = list_entry( pos2, region_t, list );
            if( r2->pattern != p2 || (r1->pivot.row == r2->pivot.row && r1->pivot.col == r2->pivot.col) )
                continue;
            
            // Compute the difference between r2's variant and variant v2
            int vv2 = (v2 - r2->variant) % base;
            vv2 = vv2 < 0 ? vv2+base : vv2;
            
            if( vv1 == vv2 &&
                r2->pivot.row - r1->pivot.row == p2_offset.row &&
                r2->pivot.col - r1->pivot.col == p2_offset.col ) {

                rfca_coord_t pivot = r1->pivot;
                // Fix the list iterators if we're removing r1 and r2
                if( tmp1 == &r2->list )
                    tmp1 = r2->list.next;
                if( tmp2 == &r1->list )
                    tmp2 = r1->list.next;
                // Remove and free both r1 and r2
                r1->pattern->usage--;
                list_del( &(r1->list) );
                region_free( r1 );
                r2->pattern->usage--;
                list_del( &(r2->list) );
                region_free( r2 );
                // Create a new region at this pivot containing p_union
                region_t* region = (region_t*)malloc( sizeof( region_t ) );
                region->pivot =pivot;
                region->pattern =p_union;
                region->variant =vv1;

                list_add_tail( &(region->list), &(v->encoded->list) );

                p_union->usage ++;
            }
        }
    }

    // If p1 and/or p2 are not used and they are not singleton patterns,
    // we remove them from the code table and free their memory.
    if( p1->usage == 0 && p1->size > 1 ) {
        list_del( &(p1->list ) );
        pattern_free( p1 );
    }
    if( p1 != p2 && p2->usage == 0 && p2->size > 1 ) {
        list_del( &(p2->list ) );
        pattern_free( p2 );
    }

    return p_union;
}

encoded_rfca_t*
encoded_create_from( rfca_t* r ) {
    // We're creating an encoded version of r using a standard code table
    encoded_rfca_t* v = (encoded_rfca_t*)malloc( sizeof( encoded_rfca_t ) );
    v->offsetCache = NULL;
    v->rfca =r;

    // The initial code table contains only one pattern
    v->codeTable = (pattern_t*)malloc( sizeof( pattern_t ) );
    INIT_LIST_HEAD( &(v->codeTable->list) );
    pattern_t* p0 = pattern_createSingle( 0 );
    list_add( &(p0->list), &(v->codeTable->list ) );

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
            region_t* region = region_create( p0, pivot );
            region->variant = value;

            // Add to the encoded dataset
            list_add_tail( &(region->list), &(v->encoded->list) );

            // Increment the pattern's usage so we can compute its code length later
            p0->usage++;
        }
    }

    // Compute the initial encoding sizes for the data and the code table
    v->stdBitsPerSingleton = -log2( 1.0 / (double)r->opts.base );
    v->stdBitsPerPivot = -log2( 1.0 / (double)v->index->nodeCount );
    v->stdBitsPerVariant = -log2( 1.0 / (double)r->opts.base );
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
    if( v->offsetCache )
        free( v->offsetCache );
    free( v );
}

static void
offsetcache_alloc( encoded_rfca_t* v ) {
    v->cacheIndex =0;
    if( !v->offsetCache ) {
        int maxOffsets = v->rfca->buffer->nodeCount / 2;
        maxOffsets *= maxOffsets;
        v->offsetCache = (pattern_offset_t*)malloc( maxOffsets * sizeof( pattern_offset_t ) );
    }
}

static bool
offsetcache_isIn( encoded_rfca_t* v, pattern_offset_t offs ) {
    for( uint64_t i =0; i < v->cacheIndex; i++ ) {
        if( offs.row == v->offsetCache[i].row &&
            offs.col == v->offsetCache[i].col )
            return true;
    }
    return false;
}

static void
offsetcache_push( encoded_rfca_t* v, pattern_offset_t offs ) {
    v->offsetCache[v->cacheIndex++] = offs;
}

int
encoded_step( encoded_rfca_t* v ) {
    // Label all the patterns so we can print them
    pattern_list_setLabels( v->codeTable ); // ONLY FOR DEBUG
    
    offsetcache_alloc( v );

    // Sort the code table by usage, then size (descending order)
    pattern_list_sortByUsageDesc( v->codeTable );

    // Take the two patterns with the largest usage
    pattern_t* p1 = list_entry( v->codeTable->list.next, pattern_t, list );
    pattern_t* p2 = list_entry( p1->list.next, pattern_t, list );

    static const int n_permutations =3;
    pattern_t* const permutations[3][2] = {
        { p1, p1 },
        { p1, p2 },
        { p2, p2 }
    };

    printf( "usages | p1: %d, p2: %d\n", p1->usage, p2->usage );

    double bestGain =0.0;
    pattern_offset_t bestP2Offset;
    int bestUsage =0, bestV1 =0, bestV2 =0, bestPerm =0;

    for( int i =0; i < n_permutations; i++ ) {
        // Prepare the patterns p1 and p2 from each permutation
        p1 = permutations[i][0];
        if( !p1 ) continue;
        p2 = permutations[i][1];
        if( !p2 ) continue;

        struct list_head *pos1, *pos2;
        list_for_each( pos1, &(v->encoded->list) ) {
            region_t* r1 = list_entry( pos1, region_t, list );
            if( r1->pattern != p1 )
                continue;
            list_for_each( pos2, &(v->encoded->list) ) { 
                region_t* r2 = list_entry( pos2, region_t, list );
                if( r2->pattern != p2 || (r1->pivot.row == r2->pivot.row && r1->pivot.col == r2->pivot.col) )
                    continue;
                
                pattern_offset_t p2_offset = pattern_offset( r1->pivot, r2->pivot );
                //if( offsetcache_isIn( v, p2_offset ) )
                //    continue;
                //offsetcache_push( v, p2_offset );

                int usage =computeUsage( v, p1, r1->variant, p2, r2->variant, p2_offset );
                double gain = computeGain( v, r1->pattern, r2->pattern, usage );
                if( gain > bestGain ) {
                    bestGain =gain;
                    bestP2Offset= p2_offset;
                    bestUsage =usage;
                    bestPerm = i;
                    bestV1 = r1->variant;
                    bestV2 = r2->variant;
                }

            }
        }
    }
    if( bestGain == 0.0 ) {
        printf( "encoded_step(): No compression gain.\n" );
        return false;
    }

    p1 = permutations[bestPerm][0];
    p2 = permutations[bestPerm][1];

    printf( "encoded_step(): merging: %c+%d and %c+%d with offset (%d,%d)\n", 
            p1->label, bestV1, p2->label, bestV2, bestP2Offset.row, bestP2Offset.col );
    printf( "encoded_step(): best usage: %d\n", bestUsage );
    printf( "encoded_step(): compression size gain: %f bits\n", bestGain );

    //updateEncoding( v, bestP );
    mergeEncodedPatterns2( v, p1, bestV1, p2, bestV2, bestP2Offset );
    pattern_list_updateCodeLength( v->codeTable, v->index->nodeCount );
    v->ctBits = computeCodeTableBits( v );
    v->encodedBits = computeEncodedBits( v );

    return true;
}

int
encoded_test( encoded_rfca_t* v ) {
    // Test piece
    
    pattern_t* p_union =NULL;
   
    region_t* r1 = list_entry( v->encoded->list.next, region_t, list );
    region_t* r2 = list_entry( r1->list.next, region_t, list );

    printf( "Code lengths | p1: %f, p2 %f; variants | p1: %d, p2: %d\n", 
            r1->pattern->codeLength, r2->pattern->codeLength, r1->variant, r2->variant );
    printf( "Total nodes: %d\n", v->index->nodeCount );

    pattern_offset_t p2_offset = pattern_offset( r1->pivot, r2->pivot );

    int usage =computeUsage( v, r1->pattern, r1->variant, r2->pattern, r2->variant, p2_offset );
    printf( "p_union usage: %d\n", usage );

  /*  int usage2 =computeUsage2( v, r1->pattern, r2->pattern, p2_offset );
    printf( "p_union usage2: %d\n", usage2 );*/

    double gain = computeGain( v, r1->pattern, r2->pattern, usage );
    printf( "Compression size gain: %f bits\n", gain );

    //updateEncoding( v, p_union );
    mergeEncodedPatterns2( v, r1->pattern, r1->variant, r2->pattern, r2->variant, p2_offset );
    pattern_list_updateCodeLength( v->codeTable, v->index->nodeCount );
    v->ctBits = computeCodeTableBits( v );
    v->encodedBits = computeEncodedBits( v );

    return 0;
}

rfca_t*
encoded_decode( encoded_rfca_t* v ) {
    rfca_t* r = rfca_create( v->rfca->opts );
    struct list_head* pos;
    list_for_each( pos, &(v->encoded->list) ) {
        region_t* region = list_entry( pos, region_t, list );

        region_apply( region, r );
    }
    return r;
}

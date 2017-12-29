/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#ifndef ENCODED_RFCA_H
#define ENCODED_RFCA_H

#include "rfca.h"
#include "region.h"
#include "pattern.h"

typedef struct {
    region_t* encoded;
    pattern_t* codeTable;
    pattern_t* singleton;
    rfca_t *rfca;
    double encodedBits;
    double ctBits;
    double stdBitsPerSingleton;
    double stdBitsPerPivot;
    double stdBitsPerVariant;
    void* offsetCache;
    uint64_t cacheIndex;
} encoded_rfca_t;

encoded_rfca_t*
encoded_createFrom( rfca_t* r );

encoded_rfca_t*
encoded_createEncodedUsing( rfca_t* r, pattern_t* codeTable );

void
encoded_free( encoded_rfca_t* v );

int
encoded_test( encoded_rfca_t* v );

int
encoded_encode( encoded_rfca_t* v );

int
encoded_encodeStep( encoded_rfca_t* v );

rfca_t*
encoded_decode( encoded_rfca_t* v );


#endif

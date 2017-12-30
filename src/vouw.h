/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#ifndef VOUW_H
#define VOUW_H

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
    double stdBitsPerOffset;
    double stdBitsPerPivot;
    double stdBitsPerVariant;
    void* buffer;
    uint64_t bufferIndex;
} vouw_t;

vouw_t*
vouw_createFrom( rfca_t* r );

vouw_t*
vouw_createEncodedUsing( rfca_t* r, pattern_t* codeTable );

void
vouw_free( vouw_t* v );

int
vouw_test( vouw_t* v );

int
vouw_encode( vouw_t* v );

int
vouw_encodeStep( vouw_t* v );

rfca_t*
vouw_decode( vouw_t* v );

#endif

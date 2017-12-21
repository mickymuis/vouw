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
    region_list_t* encoded;
    pattern_list_t* codeTable;
    rfca_buffer_t* index;
    rfca_t *rfca;
    double encodedBits;
    double ctBits;
    double stdBitsPerSingleton;
} encoded_rfca_t;

encoded_rfca_t*
encoded_create_from( rfca_t* v );

void
encoded_free( encoded_rfca_t* v );

int
encoded_step( encoded_rfca_t* v );



#endif

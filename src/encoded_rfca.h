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
    rfca_t *rfca;
    unsigned int encodedBits;
    unsigned int ctBits;
} encoded_rfca_t;

encoded_rfca_t*
encoded_create_from( rfca_t* r );

void
encoded_free( encoded_rfca_t* r );



#endif

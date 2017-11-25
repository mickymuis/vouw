/*
 * VOUW - Generation, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#ifndef RFCA_H
#define RFCA_H

#include <stdint.h>
#include <stdbool.h>

typedef struct {
    int base;
    int mode;
    uint64_t rule;
    uint8_t *input;
    int inputSize;
    int folds;
    bool right; // right-folding automaton
} rfca_opts_t;

typedef struct {
    uint8_t* cols;
    int size;
} rfca_row_t;

typedef struct {
    rfca_opts_t opts;
    int folds;
    int curRow;
    int curPos;
    uint8_t* ttable;
    rfca_row_t* rows;
    int rowCount;
    int width;
    int nodeCount;
} rfca_t;

rfca_t*
rfca_create( rfca_opts_t opts );

void
rfca_free( rfca_t* r );

uint64_t
rfca_maxRules( int mode, int base );

void
rfca_generate( rfca_t* r );




#endif


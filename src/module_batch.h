/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#ifndef MODULE_BATCH_H
#define MODULE_BATCH_H

// Defines some functions that operate on an entire rulespace

#include "rfca.h"

int module_encodeAll( rfca_opts_t opts, int argc, char** argv );

#endif



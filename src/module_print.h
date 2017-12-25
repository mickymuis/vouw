/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#ifndef MODULE_PRINT_H
#define MODULE_PRINT_H

#include "rfca.h"

void
rfca_print( rfca_t* r, bool pretty );

int module_print( rfca_opts_t opts, int argc, char** argv );

int module_printTTable( rfca_opts_t opts, int argc, char** argv );
int module_printTTable2( rfca_opts_t opts, int argc, char** argv );

#endif


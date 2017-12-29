/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#ifndef CLI_H
#define CLI_H

#include "rfca.h"

// Some boundaries and defaults
#define BASE_DEFAULT 2
#define BASE_MAX 5
#define MODE_DEFAULT 2
#define MODE_MAX 2
#define FOLDS_DEFAULT 0
#define FOLDS_MAX 10000
#define INPUT_MAX 500

void
cli_printHelp( char* exec );

bool
cli_parseOpts( rfca_opts_t *opts, char** argv[0], int* argc );

#endif


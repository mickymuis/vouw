/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include "cli.h"
#include "module.h"
#include <stdio.h>
#include <stdbool.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

void 
cli_printHelp( char* exec ) {
    fprintf( stderr, "VOUW - Generation, encoding and pattern-mining of Reduce-Fold Cellular Automata\n\
usage: %s MODULE [-r rule] [-b base][-m mode] [-i input] [--right] [-o output] [-f folds]\n\
\n\
\t MODULE        \t Specify a module to operate on the automaton (see below). \n\
\t -r \n\
\t --rule        \t Set the rule number.\n\
\t -b \n\
\t --base        \t Set the base (number of possible values per node).\n\
\t -m \n\
\t --mode        \t Sets neighborhood size.\n\
\t -i \n\
\t --input       \t Specifies the input sequence.\n\
\t -f \n\
\t --folds       \t Number of folds after reducing all input nodes.\n\
\t --right       \t Create a right-folding automaton, the default is left-folding.\n\
", exec );
    fprintf( stderr, "The following module-names are supported:\n" );
    module_printList( stderr );
}

bool
cli_parseOpts( rfca_opts_t* opts, char** argv_ptr[0], int* argc ) {

    const char* param_input = NULL;
    char **argv =*argv_ptr;

    // Parse the basic arguments from the commandline
    int i =0;
    for( i = 0; i < *argc; i++ ) {
        if( strncmp( argv[i], "-m", 2 ) == 0 || strncmp( argv[i], "--mode", 6 ) == 0 ) {
            opts->mode =atoi( argv[++i] );
        } else if( strncmp( argv[i], "-r", 2 ) == 0 || strncmp( argv[i], "--rule", 6 ) == 0 ) {
            // We use strtoull() here for 64 bit unsigned integers
            // TODO: use _strtoui64() on Windows;
            opts->rule =strtoull( argv[++i], NULL, 0 );
        } else if( strncmp( argv[i], "-b", 2 ) == 0 || strncmp( argv[i], "--base", 6 ) == 0 ) {
            opts->base =atoi( argv[++i] );
        } else if( strncmp( argv[i], "-f", 2 ) == 0 || strncmp( argv[i], "--folds", 7 ) == 0 ) {
            opts->folds =atoi( argv[++i] );
        } else if( strncmp( argv[i], "-i", 2 ) == 0 || strncmp( argv[i], "--input", 7 ) == 0 ) {
            param_input =argv[++i];
        /*} else if( strncmp( argv[i], "-o", 2 ) == 0 || strncmp( argv[i], "--output", 8 ) == 0 ) {
            param_outfile =argv[++i];*/
        } else if( strncmp( argv[i], "--right", 7 ) == 0 ) {
            opts->right =true;
        } else {
            break;
        }
    }

    // Prepare the parameters to the rfca and check the valid ranges
    if( opts->folds > FOLDS_MAX ) {
        fprintf( stderr, "Error: Parameter `folds' larger than allowed maximum (%d)\n", FOLDS_MAX );
        return false;
    }
    if( opts->mode < 2 || opts->mode > MODE_MAX ) {
        fprintf( stderr, "Error: Parameter `mode' larger than allowed maximum (%d) or smaller than 2\n", MODE_MAX );
        return false;
    }
    if( opts->base < 2 || opts->base > BASE_MAX ) {
        fprintf( stderr, "Error: Parameter `base' larger than allowed maximum (%d) or smaller than 2\n", BASE_MAX );
        return false;
    }
    uint64_t maxrules = rfca_maxRules( opts->base, opts->mode );
    if( opts->rule >= maxrules ) {
        fprintf( stderr, "Error: Parameter `rule' larger than allowed maximum (%"PRIu64") for the specified base and mode\n", maxrules );
        return false;
    }
    int inputLen = param_input == NULL ? 0 : strlen( param_input );
    if( inputLen > INPUT_MAX ) {
        fprintf( stderr, "Error: Parameter `input' longer than allowed maximum (%d)\n", INPUT_MAX );
        return false;

    }
    if( inputLen ) {
        opts->inputSize = inputLen;
        opts->input = malloc( sizeof(rfca_node_t) * opts->inputSize );
        memset( opts->input, 0, opts->inputSize );
        for( int i =0; i < inputLen; i++ ) {
            uint8_t v = param_input[i] - '0';
            if( !isdigit( param_input[i] ) || v >= opts->mode ) {
                fprintf( stderr, "Error: Parameter `input' Contains disallowed character (%c)\n", param_input[i] );
                return false;
            }
            opts->input[i] = v;
        }
    }

    *argv_ptr += i;
    *argc -= i;

    return true;
}

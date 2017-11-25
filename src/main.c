/*
 * VOUW - Generation, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include <stdio.h>
#include <stdbool.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include "vouw.h"
#include "rfca.h"

// Modules
#include "module_print.h"

// Some boundaries and defaults
#define BASE_DEFAULT 2
#define BASE_MAX 5
#define MODE_DEFAULT 2
#define MODE_MAX 2
#define FOLDS_DEFAULT 0
#define FOLDS_MAX 10000

void 
print_help( char* exec ) {
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
\t -o \n\
\t --output      \t Optional output file.\n\
\t --help        \t This help.\n\
", exec );
    fprintf( stderr, "The following module-names are supported:\n" );
    vouw_print_modules( stderr );
}

int
main( int argc, char** argv ) {

    // The basic parameters specified from the commandline
    int param_mode =MODE_DEFAULT;
    bool param_right =false;
    uint64_t param_rule =0;
    int param_base =BASE_DEFAULT;
    int param_folds =FOLD_DEFAULT;
    const char* param_outfile = NULL;
    const char* param_module = NULL;
    const char* param_input = NULL;


    // Add all module functions
    vouw_module_t module_print = {
        "print",
        "Prints the raw output generated by the automaton.",
        &print_main };
    vouw_register_module( &module_print );

    // The module is always the first argument
    if( argc == 1 ) {
        print_help( argv[0] );
        return 0;
    }
    param_module = argv[1];
    
    // Parse the basic arguments from the commandline
    for( int i = 2; i < argc; i++ ) {
        if( strncmp( argv[i], "-m", 2 ) == 0 || strncmp( argv[i], "--mode", 6 ) == 0 ) {
            param_mode =atoi( argv[++i] );
        } else if( strncmp( argv[i], "-r", 2 ) == 0 || strncmp( argv[i], "--rule", 6 ) == 0 ) {
            // We use strtoull() here for 64 bit unsigned integers
            // TODO: use _strtoui64() on Windows;
            param_rule =strtoull( argv[++i], NULL, 0 );
        } else if( strncmp( argv[i], "-b", 2 ) == 0 || strncmp( argv[i], "--base", 6 ) == 0 ) {
            param_base =atoi( argv[++i] );
        } else if( strncmp( argv[i], "-f", 2 ) == 0 || strncmp( argv[i], "--folds", 7 ) == 0 ) {
            param_folds =atoi( argv[++i] );
        } else if( strncmp( argv[i], "-i", 2 ) == 0 || strncmp( argv[i], "--input", 7 ) == 0 ) {
            param_input =argv[++i];
        } else if( strncmp( argv[i], "-o", 2 ) == 0 || strncmp( argv[i], "--output", 8 ) == 0 ) {
            param_outfile =argv[++i];
        } else if( strncmp( argv[i], "--right", 7 ) == 0 ) {
            param_right =true;
        } else if( strncmp( argv[i], "--help", 6 ) == 0 ) {
            print_help ( *argv );
            return 0;
        }
        else {
            // Just skip the rest
        }
    }

    // Prepare the parameters to the rfca and check the valid ranges
    rfca_opts_t opts;
    opts.right = param_right;
    opts.folds = param_folds;
    if( opts.folds > FOLDS_MAX ) {
        fprintf( stderr, "Error: Parameter `folds' larger than allowed maximum (%d)", FOLDS_MAX );
        return -1;
    }
    opts.folds = param_folds;
    if( opts.mode < 2 || opts.mode > MODE_MAX ) {
        fprintf( stderr, "Error: Parameter `mode' larger than allowed maximum (%d) or smaller than 2", MODE_MAX );
        return -1;
    }
    opts.base = param_base;
    if( opts.base < 2 || opts.base > BASE_MAX ) {
        fprintf( stderr, "Error: Parameter `base' larger than allowed maximum (%d) or smaller than 2", BASE_MAX );
        return -1;
    }
    opts.rule =param_rule;
    uint64_t maxrules = rfca_maxRules( opts.base, opts.mode );
    if( opts.rule >= maxrules ) {
        fprintf( stderr, "Error: Parameter `rule' larger than allowed maximum (%"PRIu64") for the specified base and mode", maxrules );
        return -1;
    }



    fprintf( stderr, "Rule: %"PRIu64", base: %d, mode %d, folds: %d\n",
            param_rule, param_base, param_mode, param_folds );

    return vouw_call_module( param_module, 0, 0, 0 );

}

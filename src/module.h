/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#ifndef MODULE_H
#define MODULE_H

#include <stdio.h>
#include "rfca.h"

// We use a simple array to hold module structures
#define MAX_MODULES 16

typedef int (*module_func_t) ( rfca_opts_t, int argc, char** argv );

typedef struct {
    const char* identifier; 
    const char* description;
    module_func_t entrypoint;
} module_t;

/* Print a list of all modules along with their descriptions */
void 
module_printList( FILE* );

/* Register a module
 */
void
module_register( module_t* m );

/* Call a module with the given identifier and arguments */
int
module_call( const char* ident, rfca_opts_t opts, int argc, char** argv );

#endif


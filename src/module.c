/*
 * VOUW - Generating, encoding and pattern-mining of Reduce-Fold Cellular Automata
 *
 * Micky Faas <micky@edukitty.org> 
 * Leiden Institute for Advanced Computer Science
 */

#include "module.h"
#include <strings.h>

static module_t* MODULES[MAX_MODULES];
static int NUM_MODULES =0;

/* Print a list of all modules along with their descriptions */
void 
module_printList( FILE* f ) {
    for( int i =0; i < NUM_MODULES; i++ ) {
        module_t* m= MODULES[i];
        fprintf( f, " `%s' \t %s\n", m->identifier, m->description );
    }
}

/* Register a module
 */
void
module_register( module_t* m ) {
    if( NUM_MODULES >= MAX_MODULES ) {
        fprintf( stderr, "Error: Number of modules exceeds MAX_MODULES\n" );
        return;
    }
    MODULES[NUM_MODULES++] = m;
}

/* Call a module with the given identifier and arguments */
int
module_call( const char* ident, rfca_opts_t opts, int argc, char** argv ) {
    module_t* m = NULL;

    for( int i =0; i < NUM_MODULES; i++ ) {
        m= MODULES[i];
        if( strcasecmp( m->identifier, ident ) == 0 )
            break;
        m =NULL;
    }
    if( m ) {
        return m->entrypoint( opts, argc, argv );
    }
    fprintf( stderr, "Error: No such module `%s'.\n", ident );
    return -1;
}


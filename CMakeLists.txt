cmake_minimum_required (VERSION 2.8)

project (vouw)

#set( CMAKE_VERBOSE_MAKEFILE ON )

list (APPEND CMAKE_C_FLAGS "-g -O2 -std=c99")

    add_executable(vouw 
        src/main.c 
        src/cli.c
        src/module.c
        src/rfca_buffer.c
        src/rfca.c
        src/ttable.c
        src/pattern.c
        src/region.c
        src/vouw.c
        src/module_print.c
        src/module_encode.c
        src/module_batch.c
	src/list_sort.c )

target_link_libraries (vouw "-lm" )

cmake_minimum_required (VERSION 2.8)

project (vouw)

#set( CMAKE_VERBOSE_MAKEFILE ON )

list (APPEND CMAKE_C_FLAGS "-g -O2 -std=c99")

    add_executable(vouw 
        src/main.c 
        src/vouw.c
        src/rfca_buffer.c
        src/rfca.c
        src/ttable.c
        src/pattern.c
        src/region.c
        src/encoded_rfca.c
        src/module_print.c
        src/module_encode.c )

target_link_libraries (vouw "-lm" )

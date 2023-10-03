#-------------------------------------------------------------------------------
# SuiteSparse/COLAMD/cmake_modules/FindCOLAMD.cmake
#-------------------------------------------------------------------------------

# The following copyright and license applies to just this file only, not to
# the library itself:
# FindCOLAMD.cmake, Copyright (c) 2022-2023, Timothy A. Davis.  All Rights Reserved.
# SPDX-License-Identifier: BSD-3-clause

#-------------------------------------------------------------------------------

# Finds the COLAMD include file and compiled library and sets:

# COLAMD_INCLUDE_DIR - where to find colamd.h
# COLAMD_LIBRARY     - dynamic COLAMD library
# COLAMD_STATIC      - static COLAMD library
# COLAMD_LIBRARIES   - libraries when using COLAMD
# COLAMD_FOUND       - true if COLAMD found

#-------------------------------------------------------------------------------

# include files for COLAMD
find_path ( COLAMD_INCLUDE_DIR
    NAMES colamd.h
    PATH_SUFFIXES include suitesparse
)

# dynamic COLAMD library (or static if no dynamic library was built)
find_library ( COLAMD_LIBRARY
    NAMES colamd colamd_static
    PATH_SUFFIXES lib build build/Release build/Debug
)

if ( MSVC )
    set ( STATIC_NAME colamd_static colamd )
else ()
    set ( STATIC_NAME colamd )
    set ( save ${CMAKE_FIND_LIBRARY_SUFFIXES} )
    set ( CMAKE_FIND_LIBRARY_SUFFIXES
        ${CMAKE_STATIC_LIBRARY_SUFFIX} ${CMAKE_FIND_LIBRARY_SUFFIXES} )
endif ()

# static COLAMD library
find_library ( COLAMD_STATIC
    NAMES ${STATIC_NAME}
    PATH_SUFFIXES lib build build/Release build/Debug
)

if ( NOT MSVC )
    # restore the CMAKE_FIND_LIBRARY_SUFFIXES variable
    set ( CMAKE_FIND_LIBRARY_SUFFIXES ${save} )
endif ()

# get version of the library from the dynamic library name
get_filename_component ( COLAMD_LIBRARY  ${COLAMD_LIBRARY} REALPATH )
get_filename_component ( COLAMD_FILENAME ${COLAMD_LIBRARY} NAME )
string (
    REGEX MATCH "[0-9]+.[0-9]+.[0-9]+"
    COLAMD_VERSION
    ${COLAMD_FILENAME}
)

# set ( COLAMD_VERSION "" )
if ( EXISTS "${COLAMD_INCLUDE_DIR}" AND NOT COLAMD_VERSION )
    # if the version does not appear in the filename, read the include file
    file ( STRINGS ${COLAMD_INCLUDE_DIR}/colamd.h COLAMD_MAJOR_STR
        REGEX "define COLAMD_MAIN_VERSION" )
    file ( STRINGS ${COLAMD_INCLUDE_DIR}/colamd.h COLAMD_MINOR_STR
        REGEX "define COLAMD_SUB_VERSION" )
    file ( STRINGS ${COLAMD_INCLUDE_DIR}/colamd.h COLAMD_PATCH_STR
        REGEX "define COLAMD_SUBSUB_VERSION" )
    #message ( STATUS "major: ${COLAMD_MAJOR_STR}" )
    #message ( STATUS "minor: ${COLAMD_MINOR_STR}" )
    #message ( STATUS "patch: ${COLAMD_PATCH_STR}" )
    string ( REGEX MATCH "[0-9]+" COLAMD_MAJOR ${COLAMD_MAJOR_STR} )
    string ( REGEX MATCH "[0-9]+" COLAMD_MINOR ${COLAMD_MINOR_STR} )
    string ( REGEX MATCH "[0-9]+" COLAMD_PATCH ${COLAMD_PATCH_STR} )
    set (COLAMD_VERSION "${COLAMD_MAJOR}.${COLAMD_MINOR}.${COLAMD_PATCH}")
endif ()

set (COLAMD_LIBRARIES ${COLAMD_LIBRARY})

include (FindPackageHandleStandardArgs)

find_package_handle_standard_args ( COLAMD
    REQUIRED_VARS COLAMD_LIBRARY COLAMD_INCLUDE_DIR
    VERSION_VAR COLAMD_VERSION
)

mark_as_advanced (
    COLAMD_INCLUDE_DIR
    COLAMD_LIBRARY
    COLAMD_STATIC
    COLAMD_LIBRARIES
)

if ( COLAMD_FOUND )
  if ( NOT TARGET SuiteSparse::COLAMD )
    add_library(SuiteSparse::COLAMD INTERFACE IMPORTED GLOBAL)
    set_target_properties(SuiteSparse::COLAMD PROPERTIES INTERFACE_LINK_LIBRARIES "${COLAMD_LIBRARIES}")
  endif ()
else ()
    message ( STATUS "COLAMD not found" )
    set ( COLAMD_INCLUDE_DIR "" )
    set ( COLAMD_LIBRARIES "" )
    set ( COLAMD_LIBRARY "" )
    set ( COLAMD_STATIC "" )
endif ()


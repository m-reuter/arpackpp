#-------------------------------------------------------------------------------
# SuiteSparse/SuiteSparse_config/cmake_modules/FindSuiteSparse_config.cmake
#-------------------------------------------------------------------------------

# The following copyright and license applies to just this file only, not to
# the library itself:
# FindSuiteSparse_config.cmake, Copyright (c) 2022-2023, Timothy A. Davis.  All Rights Reserved.
# SPDX-License-Identifier: BSD-3-clause

#-------------------------------------------------------------------------------

# Finds the SuiteSparse_config include file and compiled library and sets:

# SUITESPARSE_CONFIG_INCLUDE_DIR - where to find SuiteSparse_config.h
# SUITESPARSE_CONFIG_LIBRARY     - dynamic SuiteSparse_config library
# SUITESPARSE_CONFIG_STATIC      - static SuiteSparse_config library
# SUITESPARSE_CONFIG_LIBRARIES   - libraries when using SuiteSparse_config
# SUITESPARSE_CONFIG_FOUND       - true if SuiteSparse_config found

#-------------------------------------------------------------------------------

# include files for SuiteSparse_config
find_path ( SUITESPARSE_CONFIG_INCLUDE_DIR
    NAMES SuiteSparse_config.h
    PATH_SUFFIXES include suitesparse
)

# dynamic SuiteSparse_config (or static if no dynamic library was built)
find_library ( SUITESPARSE_CONFIG_LIBRARY
    NAMES suitesparseconfig suitesparseconfig_static
    PATH_SUFFIXES lib build build/Release build/Debug
)

if ( MSVC )
    set ( STATIC_NAME suitesparseconfig_static suitesparseconfig )
else ()
    set ( STATIC_NAME suitesparseconfig )
    set ( save ${CMAKE_FIND_LIBRARY_SUFFIXES} )
    message ( STATUS "original library suffixes: ${CMAKE_FIND_LIBRARY_SUFFIXES}" )
    set ( CMAKE_FIND_LIBRARY_SUFFIXES
        ${CMAKE_STATIC_LIBRARY_SUFFIX} ${CMAKE_FIND_LIBRARY_SUFFIXES} )
    message ( STATUS "revised for static search: ${CMAKE_FIND_LIBRARY_SUFFIXES}" )
endif ()

# static libraries for SuiteSparse_config
find_library ( SUITESPARSE_CONFIG_STATIC
    NAMES ${STATIC_NAME}
    PATH_SUFFIXES lib build build/Release build/Debug
)

if ( NOT MSVC )
    # restore the CMAKE_FIND_LIBRARY_SUFFIXES variable
    set ( CMAKE_FIND_LIBRARY_SUFFIXES ${save} )
endif ()

# get version of the library from the dynamic library filename, if present
get_filename_component ( SUITESPARSE_CONFIG_LIBRARY  ${SUITESPARSE_CONFIG_LIBRARY} REALPATH )
get_filename_component ( SUITESPARSE_CONFIG_FILENAME ${SUITESPARSE_CONFIG_LIBRARY} NAME )
string (
    REGEX MATCH "[0-9]+.[0-9]+.[0-9]+"
    SUITESPARSE_CONFIG_VERSION
    ${SUITESPARSE_CONFIG_FILENAME}
)

# set ( SUITESPARSE_CONFIG_VERSION "" )
if ( EXISTS "${SUITESPARSE_CONFIG_INCLUDE_DIR}" AND NOT SUITESPARSE_CONFIG_VERSION )
    # if the version does not appear in the filename, read the include file
    file ( STRINGS ${SUITESPARSE_CONFIG_INCLUDE_DIR}/SuiteSparse_config.h SUITESPARSE_CONFIG_MAJOR_STR
        REGEX "define SUITESPARSE_MAIN_VERSION" )
    file ( STRINGS ${SUITESPARSE_CONFIG_INCLUDE_DIR}/SuiteSparse_config.h SUITESPARSE_CONFIG_MINOR_STR
        REGEX "define SUITESPARSE_SUB_VERSION" )
    file ( STRINGS ${SUITESPARSE_CONFIG_INCLUDE_DIR}/SuiteSparse_config.h SUITESPARSE_CONFIG_PATCH_STR
        REGEX "define SUITESPARSE_SUBSUB_VERSION" )
    #message ( STATUS "major: ${SUITESPARSE_CONFIG_MAJOR_STR}" )
    #message ( STATUS "minor: ${SUITESPARSE_CONFIG_MINOR_STR}" )
    #message ( STATUS "patch: ${SUITESPARSE_CONFIG_PATCH_STR}" )
    string ( REGEX MATCH "[0-9]+" SUITESPARSE_CONFIG_MAJOR ${SUITESPARSE_CONFIG_MAJOR_STR} )
    string ( REGEX MATCH "[0-9]+" SUITESPARSE_CONFIG_MINOR ${SUITESPARSE_CONFIG_MINOR_STR} )
    string ( REGEX MATCH "[0-9]+" SUITESPARSE_CONFIG_PATCH ${SUITESPARSE_CONFIG_PATCH_STR} )
    set (SUITESPARSE_CONFIG_VERSION "${SUITESPARSE_CONFIG_MAJOR}.${SUITESPARSE_CONFIG_MINOR}.${SUITESPARSE_CONFIG_PATCH}")
endif ()

# libaries when using SuiteSparse_config
set (SUITESPARSE_CONFIG_LIBRARIES ${SUITESPARSE_CONFIG_LIBRARY})

include ( FindPackageHandleStandardArgs )

find_package_handle_standard_args ( SuiteSparse_config
    REQUIRED_VARS SUITESPARSE_CONFIG_LIBRARY SUITESPARSE_CONFIG_INCLUDE_DIR
    VERSION_VAR SUITESPARSE_CONFIG_VERSION
    REASON_FAILURE_MESSAGE result
)

mark_as_advanced (
    SUITESPARSE_CONFIG_INCLUDE_DIR
    SUITESPARSE_CONFIG_LIBRARY
    SUITESPARSE_CONFIG_STATIC
    SUITESPARSE_CONFIG_LIBRARIES
)

if ( SUITESPARSE_CONFIG_FOUND )
  if ( NOT TARGET SuiteSparse::SuiteSparse_config )
    add_library(SuiteSparse::SuiteSparse_config INTERFACE IMPORTED GLOBAL)
    set_target_properties(SuiteSparse::SuiteSparse_config PROPERTIES INTERFACE_LINK_LIBRARIES "${SUITESPARSE_CONFIG_LIBRARIES}")
  endif ()
else ()
    message ( STATUS "SuiteSparse_config not found" )
    set ( SUITESPARSE_CONFIG_INCLUDE_DIR "" )
    set ( SUITESPARSE_CONFIG_LIBRARIES "" )
    set ( SUITESPARSE_CONFIG_LIBRARY "" )
    set ( SUITESPARSE_CONFIG_STATIC "" )
endif ()



#-------------------------------------------------------------------------------
# SuiteSparse/CAMD/cmake_modules/FindCAMD.cmake
#-------------------------------------------------------------------------------

# The following copyright and license applies to just this file only, not to
# the library itself:
# FindCAMD.cmake, Copyright (c) 2022-2023, Timothy A. Davis.  All Rights Reserved.
# SPDX-License-Identifier: BSD-3-clause

#-------------------------------------------------------------------------------

# Finds the CAMD include file and compiled library and sets:

# CAMD_INCLUDE_DIR - where to find camd.h
# CAMD_LIBRARY     - dynamic CAMD library
# CAMD_STATIC      - static CAMD library
# CAMD_LIBRARIES   - libraries when using CAMD
# CAMD_FOUND       - true if CAMD found

#-------------------------------------------------------------------------------

# include files for CAMD
find_path ( CAMD_INCLUDE_DIR
    NAMES camd.h
    PATH_SUFFIXES include suitesparse
)

# dynamic CAMD library (or static if no dynamic library was built)
find_library ( CAMD_LIBRARY
    NAMES camd camd_static
    PATH_SUFFIXES lib build build/Release build/Debug
)

if ( MSVC )
    set ( STATIC_NAME camd_static camd )
else ()
    set ( STATIC_NAME camd )
    set ( save ${CMAKE_FIND_LIBRARY_SUFFIXES} )
    set ( CMAKE_FIND_LIBRARY_SUFFIXES
        ${CMAKE_STATIC_LIBRARY_SUFFIX} ${CMAKE_FIND_LIBRARY_SUFFIXES} )
endif ()

# static CAMD library
find_library ( CAMD_STATIC
    NAMES ${STATIC_NAME}
    PATH_SUFFIXES lib build build/Release build/Debug
)

if ( NOT MSVC )
    # restore the CMAKE_FIND_LIBRARY_SUFFIXES variable
    set ( CMAKE_FIND_LIBRARY_SUFFIXES ${save} )
endif ()

# get version of the library from the dynamic library filename
get_filename_component ( CAMD_LIBRARY  ${CAMD_LIBRARY} REALPATH )
get_filename_component ( CAMD_FILENAME ${CAMD_LIBRARY} NAME )
string (
    REGEX MATCH "[0-9]+.[0-9]+.[0-9]+"
    CAMD_VERSION
    ${CAMD_FILENAME}
)

# set ( CAMD_VERSION "" )
if ( EXISTS "${CAMD_INCLUDE_DIR}" AND NOT CAMD_VERSION )
    # if the version does not appear in the filename, read the include file
    file ( STRINGS ${CAMD_INCLUDE_DIR}/camd.h CAMD_MAJOR_STR
        REGEX "define CAMD_MAIN_VERSION" )
    file ( STRINGS ${CAMD_INCLUDE_DIR}/camd.h CAMD_MINOR_STR
        REGEX "define CAMD_SUB_VERSION" )
    file ( STRINGS ${CAMD_INCLUDE_DIR}/camd.h CAMD_PATCH_STR
        REGEX "define CAMD_SUBSUB_VERSION" )
    #message ( STATUS "major: ${CAMD_MAJOR_STR}" )
    #message ( STATUS "minor: ${CAMD_MINOR_STR}" )
    #message ( STATUS "patch: ${CAMD_PATCH_STR}" )
    string ( REGEX MATCH "[0-9]+" CAMD_MAJOR ${CAMD_MAJOR_STR} )
    string ( REGEX MATCH "[0-9]+" CAMD_MINOR ${CAMD_MINOR_STR} )
    string ( REGEX MATCH "[0-9]+" CAMD_PATCH ${CAMD_PATCH_STR} )
    set (CAMD_VERSION "${CAMD_MAJOR}.${CAMD_MINOR}.${CAMD_PATCH}")
endif ()

set ( CAMD_LIBRARIES ${CAMD_LIBRARY} )

include (FindPackageHandleStandardArgs)

find_package_handle_standard_args ( CAMD
    REQUIRED_VARS CAMD_LIBRARY CAMD_INCLUDE_DIR
    VERSION_VAR CAMD_VERSION
)

mark_as_advanced (
    CAMD_INCLUDE_DIR
    CAMD_LIBRARY
    CAMD_STATIC
    CAMD_LIBRARIES
)

if ( CAMD_FOUND )
  if ( NOT TARGET SuiteSparse::CAMD )
    add_library(SuiteSparse::CAMD INTERFACE IMPORTED GLOBAL)
    set_target_properties(SuiteSparse::CAMD PROPERTIES INTERFACE_LINK_LIBRARIES "${CAMD_LIBRARIES}")
  endif ()
else ()
    message ( STATUS "CAMD not found" )
    set ( CAMD_INCLUDE_DIR "" )
    set ( CAMD_LIBRARIES "" )
    set ( CAMD_LIBRARY "" )
    set ( CAMD_STATIC "" )
endif ()


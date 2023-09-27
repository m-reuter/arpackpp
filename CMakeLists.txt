cmake_minimum_required (VERSION 3.18)

project(arpackpp VERSION 2.3.0
                 DESCRIPTION "ARPACK++: C++ interface to ARPACK"
                 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(arpackpp_MAJOR_VERSION 2)
set(arpackpp_MINOR_VERSION 3)
set(arpackpp_PATCH_VERSION 0)
set(arpackpp_VERSION ${arpackpp_MAJOR_VERSION}.${arpackpp_MINOR_VERSION}.${arpackpp_PATCH_VERSION})

option(ENABLE_EXAMPLES "Enable examples" ON)
option(ENABLE_FORTRAN "Enable Fortran language (for static linking of ARPACK)" OFF)
option(ENABLE_SUPERLU "Enable SUPERLU" OFF)
option(ENABLE_UMFPACK "Enable UMFPACK" OFF)
option(ENABLE_CHOLMOD "Enable CHOLMOD" OFF)

set (suitesparse_static 0)

if (PREFER_SUITESPARSE_STATIC)
  set (suitesparse_static 1)
endif()

if (ENABLE_FORTRAN)
  enable_language(Fortran)
endif()

# In case the dependencies are installed locally using the shell scripts, we need to
# let CMake know, where to look for the config-file packages.
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${CMAKE_SOURCE_DIR}/external)

# In case dependencies were installed by a package manager, we need to let CMake know
# where to look for the modules.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules)

# Custom function to find packages. First try config mode. If the target isn't found,
# use module mode.
function(find_package_custom name target)
  find_package(${name} CONFIG QUIET)
  
  if (NOT TARGET ${target})
    find_package(${name} REQUIRED)
  else()
    get_target_property(location ${target} LOCATION)
    message( STATUS "Found ${name} [config]: ${location}" )
  endif()
endfunction()

# Find BLAS and LAPACK

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Find arpack-ng

find_package_custom(arpackng ARPACK::ARPACK)

get_target_property(type ARPACK::ARPACK TYPE)
if (NOT ENABLE_FORTRAN AND NOT WIN32 AND type MATCHES "^STATIC")
  message( STATUS "Enabling Fortran for static linking to ARPACK..." )
  enable_language(Fortran)
endif()

# Find SuperLU

if (ENABLE_SUPERLU)
  find_package_custom(superlu superlu::superlu)
endif()

# Find SuiteSparse

if (ENABLE_CHOLMOD OR ENABLE_UMFPACK)
  if (${suitesparse_static})
    message(STATUS "Linking SuiteSparse static")
    find_package_custom ( SuiteSparse_config SuiteSparse::SuiteSparse_config_static )
    find_package_custom ( AMD SuiteSparse::AMD_static )
    find_package_custom ( CAMD SuiteSparse::CAMD_static )
    find_package_custom ( CCOLAMD SuiteSparse::CCOLAMD_static )
    find_package_custom ( COLAMD SuiteSparse::COLAMD_static )
    find_package_custom ( CHOLMOD SuiteSparse::CHOLMOD_static )

    if (NOT TARGET SuiteSparse::CHOLMOD_static)
      message(FATAL_ERROR "CHOLMOD_static not found")
    endif()
  else()
    find_package_custom ( SuiteSparse_config SuiteSparse::SuiteSparse_config )
    find_package_custom ( AMD SuiteSparse::AMD )
    find_package_custom ( CAMD SuiteSparse::CAMD )
    find_package_custom ( CCOLAMD SuiteSparse::CCOLAMD )
    find_package_custom ( COLAMD SuiteSparse::COLAMD )
    find_package_custom ( CHOLMOD SuiteSparse::CHOLMOD )

    if (NOT TARGET SuiteSparse::CHOLMOD)
      message(FATAL_ERROR "CHOLMOD not found")
    endif()
  endif()
endif()

if (ENABLE_UMFPACK)
  if (${suitesparse_static})
    find_package_custom ( UMFPACK SuiteSparse::UMFPACK_static )

    if (NOT TARGET SuiteSparse::UMFPACK_static )
      message(FATAL_ERROR "UMFPACK_static not found")
    endif()
  else()
    find_package_custom ( UMFPACK SuiteSparse::UMFPACK )

    if (NOT TARGET SuiteSparse::UMFPACK )
      message(FATAL_ERROR "UMFPACK not found")
    endif()
  endif()
endif()

# ARPACK++ target

add_library(arpackpp INTERFACE)
target_link_libraries(arpackpp INTERFACE ARPACK::ARPACK)
target_include_directories(arpackpp INTERFACE include)

# Examples

if (ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

# Install

file (GLOB HEADERS "include/*.h")
install (FILES ${HEADERS} DESTINATION include/arpackpp)

file (GLOB_RECURSE EXAMPLE_HEADERS "examples/*.h")
install (FILES ${EXAMPLE_HEADERS} DESTINATION include/arpackpp)

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
option(ENABLE_SUPERLU "Enable SUPERLU" OFF)
option(ENABLE_UMFPACK "Enable UMFPACK" OFF)
option(ENABLE_CHOLMOD "Enable CHOLMOD" OFF)

# In case the dependencies are installed with the --local-install option, we need to
# let CMake know, where to look for the config-file packages.
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${CMAKE_SOURCE_DIR}/external)

# Find BLAS

find_package(BLAS REQUIRED)

if (NOT BLAS_FOUND)
  message(FATAL_ERROR "BLAS not found")
endif()

# Find LAPACK

find_package(LAPACK REQUIRED)

if (NOT LAPACK_FOUND)
  message(FATAL_ERROR "LAPACK not found")
endif()

# Find arpack-ng

find_package (arpackng REQUIRED)

if (TARGET ARPACK::ARPACK)
  get_target_property(location ARPACK::ARPACK LOCATION)
  message( STATUS "Found ARPACK: ${location} (version ${arpackng_VERSION})" )
else()
  message(FATAL_ERROR "ARPACK not found")
endif()

# Find SuperLU

if (ENABLE_SUPERLU)
  find_package(superlu REQUIRED)

  if (TARGET superlu::superlu)
    get_target_property(location superlu::superlu LOCATION)
    message( STATUS "Found SuperLU: ${location} (version ${superlu_VERSION})" )
  else()
    message(FATAL_ERROR "SuperLU not found")
  endif()
endif()

# Find SuiteSparse

if (ENABLE_CHOLMOD OR ENABLE_UMFPACK)
  find_package ( SuiteSparse_config REQUIRED )
  find_package ( AMD REQUIRED )
  find_package ( CAMD REQUIRED )
  find_package ( CCOLAMD REQUIRED )
  find_package ( COLAMD REQUIRED )
  find_package ( CHOLMOD REQUIRED )

  if (NOT TARGET SuiteSparse::CHOLMOD)
    message(FATAL_ERROR "CHOLMOD not found")
  endif()
endif()

if (ENABLE_UMFPACK)
  find_package (UMFPACK REQUIRED)

  if (NOT TARGET SuiteSparse::UMFPACK)
    message(FATAL_ERROR "UMFPACK not found")
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

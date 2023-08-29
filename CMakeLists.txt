cmake_minimum_required (VERSION 3.16)

project(arpackpp VERSION 2.3.0
                 DESCRIPTION "ARPACK++: C++ interface to ARPACK"
                 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(arpackpp_MAJOR_VERSION 2)
set(arpackpp_MINOR_VERSION 3)
set(arpackpp_PATCH_VERSION 0)
set(arpackpp_VERSION ${arpackpp_MAJOR_VERSION}.${arpackpp_MINOR_VERSION}.${arpackpp_PATCH_VERSION})

option(ENABLE_SUPERLU "Enable SUPERLU" OFF)
option(ENABLE_UMFPACK "Enable UMFPACK" OFF)
option(ENABLE_CHOLMOD "Enable CHOLMOD" OFF)
option(ENABLE_EXAMPLES "Enable examples" ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/SuiteSparse/")

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

# ARPACK

find_package (ARPACK)

if (NOT ARPACK_FOUND)
  message(FATAL_ERROR "ARPACK not found")
endif()

# SUPERLU

if (ENABLE_SUPERLU)
  find_package(superlu REQUIRED)

  if (TARGET superlu::superlu)
    message( STATUS "Found SuperLU: version ${superlu_VERSION}" )
  else()
    message(FATAL_ERROR "SuperLU not found")
  endif()
endif()

# Suitesparse

set(SUITESPARSE_INCLUDE_DIRS "")
set(SUITESPARSE_LIBRARIES "")

function(suitesparse_load name)
  find_package(${name} REQUIRED)

  string(TOUPPER "${name}" name_UPPER)

  if (NOT ${${name_UPPER}_FOUND})
    message(FATAL_ERROR "${name} not found")
  endif()
  
  if (NOT ${${name_UPPER}_INCLUDE_DIR} IN_LIST SUITESPARSE_INCLUDE_DIRS)
    set(SUITESPARSE_INCLUDE_DIRS ${SUITESPARSE_INCLUDE_DIRS} ${${name_UPPER}_INCLUDE_DIR} PARENT_SCOPE)
  endif()
  
  set(SUITESPARSE_LIBRARIES ${SUITESPARSE_LIBRARIES} ${${name_UPPER}_LIBRARY} PARENT_SCOPE)
endfunction()

if (ENABLE_CHOLMOD OR ENABLE_UMFPACK)
  suitesparse_load(SuiteSparse_config)
  suitesparse_load(AMD)
  suitesparse_load(CAMD)
  suitesparse_load(COLAMD)
  suitesparse_load(CCOLAMD)
  suitesparse_load(CHOLMOD)

  if (ENABLE_UMFPACK)
    suitesparse_load(UMFPACK)
  endif()
endif()

# ARPACK++ target

add_library(arpackpp INTERFACE)
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

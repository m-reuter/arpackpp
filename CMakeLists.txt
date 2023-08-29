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

# Find BLAS

find_package(BLAS REQUIRED)

if (NOT TARGET BLAS::BLAS) # Create target if cmake version < 3.18
  add_library(BLAS::BLAS INTERFACE IMPORTED)
  set_target_properties(BLAS::BLAS PROPERTIES INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}")
endif()

if (BLAS_LIBRARIES)
  message( STATUS "BLAS found: ${BLAS_LIBRARIES}" )
endif()

# Find LAPACK

find_package(LAPACK REQUIRED)

if (NOT TARGET LAPACK::LAPACK) # Create target if cmake version < 3.18
  add_library(LAPACK::LAPACK INTERFACE IMPORTED)
  set_target_properties(LAPACK::LAPACK PROPERTIES INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}")
endif()

if (LAPACK_LIBRARIES)
  message( STATUS "LAPACK found: ${LAPACK_LIBRARIES}" )
endif()

# ARPACK

find_library (ARPACK_LIB NAMES arpack PATHS external external/ARPACK)

if (NOT ARPACK_LIB)
  message(SEND_ERROR "libarpack is required but could not be found")
  set(ABORT_CONFIG TRUE)
endif()
message( STATUS "ARPACK found: ${ARPACK_LIB}" )

# SUPERLU

if (ENABLE_SUPERLU)
  find_path(SUPERLU_INC supermatrix.h "external/SuperLU/SRC" "/usr/local/include/superlu" )
  if (NOT SUPERLU_INC)
    message(STATUS "SuperLU (supermatrix.h) is required but could not be found")
    set(ABORT_CONFIG TRUE)
  endif()
  
  find_library(SUPERLU_LIB superlu.a libsuperlu.a  
    PATHS ${SUPERLU_INC} ${SUPERLU_INC}/../ external /usr/local/lib )
  if (NOT SUPERLU_LIB)
    message(STATUS "SuperLU Library is required but could not be found")
    set(ABORT_CONFIG TRUE)
  endif()
  include_directories( ${SUPERLU_INC} )
  message( STATUS "SUPERLU_LIB found: ${SUPERLU_LIB}" ) 
endif()

# Suitesparse

if (ENABLE_UMFPACK OR ENABLE_CHOLMOD)
  # Suitesparse CHOLMOD
  find_path(SUITESPARSE_DIR SuiteSparse_config/SuiteSparse_config.h "external/SuiteSparse")
  if (NOT SUITESPARSE_DIR)
    # nothing in external, search globally for libraries and headers
    find_library(CHOLMOD_LIB libcholmod.a)
    find_library(COLAMD_LIB libcolamd.a)
    find_library(CCOLAMD_LIB libccolamd.a)
    find_library(METIS_LIB libmetis.a)
    if (NOT METIS_LIB)
      find_library(METIS_LIB libmetis.dylib)
    endif()
    if (NOT METIS_LIB)
      find_library(METIS_LIB libmetis.so)
    endif()
    find_library(CAMD_LIB libcamd.a)
    find_library(AMD_LIB libamd.a)

    # Suitesparse UMFPACK
    find_library(UMFPACK_LIB libumfpack.a)
    find_library(SUITESPARSE_LIB libsuitesparseconfig.a)
   
    find_path( CHOLMOD_INCLUDE_DIR cholmod.h
               PATHS /usr/local/include
               /usr/include
               /usr/include/suitesparse/
               ${CMAKE_SOURCE_DIR}/MacOS/Libs/cholmod
               PATH_SUFFIXES cholmod/ CHOLMOD/ )
    if (NOT CHOLMOD_INCLUDE_DIR)
      message(STATUS "SuiteSparse chomod.h is required but could not be found")
      set(ABORT_CONFIG TRUE)
    endif()

    include_directories(
      ${CHOLMOD_INCLUDE_DIR}
    )
    
    
  else (NOT SUITESPARSE_DIR)
    # found suite sparse locally in external dir
    message( STATUS "SUITESPARSE_DIR found: ${SUITESPARSE_DIR}" )

    find_library(CHOLMOD_LIB libcholmod.a ${SUITESPARSE_DIR}/CHOLMOD/Lib  )
    find_library(COLAMD_LIB libcolamd.a ${SUITESPARSE_DIR}/COLAMD/Lib  )
    find_library(CCOLAMD_LIB libccolamd.a ${SUITESPARSE_DIR}/CCOLAMD/Lib  )
    find_library(METIS_LIB libmetis.a ${SUITESPARSE_DIR}/metis-4.0 )
    if (NOT METIS_LIB)
      find_library(METIS_LIB libmetis.dylib)
    endif()
    if (NOT METIS_LIB)
      find_library(METIS_LIB libmetis.so)
    endif()
    find_library(CAMD_LIB libcamd.a ${SUITESPARSE_DIR}/CAMD/Lib )
    find_library(AMD_LIB libamd.a ${SUITESPARSE_DIR}/AMD/Lib )

    # Suitesparse UMFPACK
    find_library(UMFPACK_LIB libumfpack.a ${SUITESPARSE_DIR}/UMFPACK/Lib  )
    find_library(SUITESPARSE_LIB libsuitesparseconfig.a ${SUITESPARSE_DIR}/SuiteSparse_config )

    include_directories(
      ${SUITESPARSE_DIR}/CHOLMOD/Include
      ${SUITESPARSE_DIR}/SuiteSparse_config
    )


  endif()
  
  set(SSPARSE_LIBS
    ${UMFPACK_LIB}
    ${CHOLMOD_LIB}
    ${COLAMD_LIB}
    ${CCOLAMD_LIB}
    ${METIS_LIB}
    ${CAMD_LIB}
    ${AMD_LIB}
    ${SUITESPARSE_LIB} 
  )
  
  # some linux systems need specifically librt (for SuiteSparse):
  if (NOT APPLE)
    set (SSPARSE_LIBS ${SSPARSE_LIBS} -lrt)
  endif()
  
endif()

# ARPACK++ target

add_custom_target( Message ALL
  COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "Run either \\\"make examples\\\" or \\\"make install\\\"."
  COMMENT "Header Library: Nothing to compile!" )

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

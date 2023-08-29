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

option(SUPERLU "Enable SUPERLU" OFF)
option(UMFPACK "Enable UMFPACK" OFF)
option(CHOLMOD "Enable CHOLMOD" OFF)

add_custom_target (examples)

function(examples list_name)
  foreach(l ${${list_name}})
    get_filename_component(lwe ${l} NAME_WE)
    add_executable(${lwe} EXCLUDE_FROM_ALL ${l})
    target_link_libraries(${lwe}
      PUBLIC
        ${ARPACK_LIB}
        $<BUILD_INTERFACE:LAPACK::LAPACK>
        $<BUILD_INTERFACE:BLAS::BLAS>)
    add_dependencies (examples ${lwe})
  endforeach()
endfunction(examples)

function(examplesslu list_name)
  foreach(l ${${list_name}})
    get_filename_component(lwe ${l} NAME_WE)
    add_executable(${lwe} EXCLUDE_FROM_ALL ${l})
    target_link_libraries(${lwe}
      PUBLIC
        ${ARPACK_LIB}
        ${SUPERLU_LIB}
        $<BUILD_INTERFACE:LAPACK::LAPACK>
        $<BUILD_INTERFACE:BLAS::BLAS>)
    add_dependencies (examples ${lwe})
  endforeach()
endfunction(examplesslu)

function(examplesls list_name)
  foreach(l ${${list_name}})
    get_filename_component(lwe ${l} NAME_WE)
    add_executable(${lwe} EXCLUDE_FROM_ALL ${l})
    target_link_libraries(${lwe}
      PUBLIC
        ${ARPACK_LIB}
        ${SSPARSE_LIBS}
        $<BUILD_INTERFACE:LAPACK::LAPACK>
        $<BUILD_INTERFACE:BLAS::BLAS>)
    add_dependencies (examples ${lwe})
  endforeach()
endfunction(examplesls)

function(make_global_path invar cstr)
  set (infile ${${invar}})
  if (infile)
   if (EXISTS ${infile})
    if(NOT IS_ABSOLUTE ${infile})
      get_filename_component(_infile_tmp ${CMAKE_BINARY_DIR}/${infile} REALPATH)
      set(${invar}  ${_infile_tmp} CACHE STRING ${cstr} FORCE)
    endif()
  else()
    MESSAGE(STATUS "${infile} does not exist!")
    set(${invar} ${invar}-NOTFOUND CACHE STRING ${cstr} FORCE)   
  endif()
 endif()
endfunction(make_global_path)


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


IF (NOT ARPACK_LIB)
  MESSAGE(SEND_ERROR "libarpack is required but could not be found")
  SET(ABORT_CONFIG TRUE)
ENDIF()
make_global_path(ARPACK_LIB "ARPACK Fortran Library")
message( STATUS "ARPACK_LIB found: ${ARPACK_LIB}" )

# SUPERLU
if (SUPERLU)
  find_path(SUPERLU_INC supermatrix.h "external/SuperLU/SRC" "/usr/local/include/superlu" )
  IF (NOT SUPERLU_INC)
    MESSAGE(STATUS "SuperLU (supermatrix.h) is required but could not be found")
    SET(ABORT_CONFIG TRUE)
  ENDIF()
  make_global_path(SUPERLU_INC "SUPERLU Include Path")
  
  find_library(SUPERLU_LIB superlu.a libsuperlu.a  
    PATHS ${SUPERLU_INC} ${SUPERLU_INC}/../ external /usr/local/lib )
  IF (NOT SUPERLU_LIB)
    MESSAGE(STATUS "SuperLU Library is required but could not be found")
    SET(ABORT_CONFIG TRUE)
  ENDIF()
  make_global_path(SUPERLU_LIB "SUPERLU Library")
  include_directories( ${SUPERLU_INC} )
  message( STATUS "SUPERLU_LIB found: ${SUPERLU_LIB}" ) 
endif()

# Suitesparse 
if (UMFPACK OR CHOLMOD)
  # Suitesparse CHOLMOD
  find_path(SUITESPARSE_DIR SuiteSparse_config/SuiteSparse_config.h "external/SuiteSparse")
  IF (NOT SUITESPARSE_DIR)
    # nothing in external, search globally for libraries and headers
    find_library(CHOLMOD_LIB libcholmod.a)
    find_library(COLAMD_LIB libcolamd.a)
    find_library(CCOLAMD_LIB libccolamd.a)
    find_library(METIS_LIB libmetis.a)
    IF (NOT METIS_LIB)
      find_library(METIS_LIB libmetis.dylib)
    ENDIF()
    IF (NOT METIS_LIB)
      find_library(METIS_LIB libmetis.so)
    ENDIF()
    find_library(CAMD_LIB libcamd.a)
    find_library(AMD_LIB libamd.a)

    # Suitesparse UMFPACK
    find_library(UMFPACK_LIB libumfpack.a)
    find_library(SUITESPARSE_LIB libsuitesparseconfig.a)
   
    FIND_PATH( CHOLMOD_INCLUDE_DIR cholmod.h
        	      PATHS /usr/local/include 
        	      /usr/include 
        	      /usr/include/suitesparse/ 
        	      ${CMAKE_SOURCE_DIR}/MacOS/Libs/cholmod
              	PATH_SUFFIXES cholmod/ CHOLMOD/ )
    IF (NOT CHOLMOD_INCLUDE_DIR)
      MESSAGE(STATUS "SuiteSparse chomod.h is required but could not be found")
      SET(ABORT_CONFIG TRUE)
    ENDIF()

    include_directories(
      ${CHOLMOD_INCLUDE_DIR}
    )
    
    
  ELSE (NOT SUITESPARSE_DIR)
    # found suite sparse locally in external dir
    make_global_path(SUITESPARSE_DIR "SuiteSparse Directory")
    message( STATUS "SUITESPARSE_DIR found: ${SUITESPARSE_DIR}" )

    find_library(CHOLMOD_LIB libcholmod.a ${SUITESPARSE_DIR}/CHOLMOD/Lib  )
    find_library(COLAMD_LIB libcolamd.a ${SUITESPARSE_DIR}/COLAMD/Lib  )
    find_library(CCOLAMD_LIB libccolamd.a ${SUITESPARSE_DIR}/CCOLAMD/Lib  )
    find_library(METIS_LIB libmetis.a ${SUITESPARSE_DIR}/metis-4.0 )
    IF (NOT METIS_LIB)
      find_library(METIS_LIB libmetis.dylib)
    ENDIF()
    IF (NOT METIS_LIB)
      find_library(METIS_LIB libmetis.so)
    ENDIF()
    find_library(CAMD_LIB libcamd.a ${SUITESPARSE_DIR}/CAMD/Lib )
    find_library(AMD_LIB libamd.a ${SUITESPARSE_DIR}/AMD/Lib )

    # Suitesparse UMFPACK
    find_library(UMFPACK_LIB libumfpack.a ${SUITESPARSE_DIR}/UMFPACK/Lib  )
    find_library(SUITESPARSE_LIB libsuitesparseconfig.a ${SUITESPARSE_DIR}/SuiteSparse_config )

    include_directories(
      ${SUITESPARSE_DIR}/CHOLMOD/Include
      ${SUITESPARSE_DIR}/SuiteSparse_config
    )


  ENDIF()
  
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

include_directories(
  include
  examples/matrices/complex
  examples/matrices/nonsym
  examples/matrices/sym
  examples/matprod
  examples/matprod/complex
  examples/matprod/nonsym
  examples/matprod/sym
  examples/areig
  
)

file (GLOB HEADERS "include/*.h")
install (FILES ${HEADERS} DESTINATION include/arpackpp)

file (GLOB_RECURSE EXAMPLE_HEADERS "examples/*.h")
install (FILES ${EXAMPLE_HEADERS} DESTINATION include/arpackpp)

# examples product
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/product/)
file(GLOB product_complex examples/product/complex/*.cc)
examples(product_complex)
file(GLOB product_nonsym examples/product/nonsym/*.cc)
examples(product_nonsym)
file(GLOB product_simple examples/product/simple/*.cc)
examples(product_simple)
file(GLOB product_sym examples/product/sym/*.cc)
examples(product_sym)

# examples reverse
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/reverse/)
file(GLOB reverse_complex examples/reverse/complex/*.cc)
examples(reverse_complex)
file(GLOB reverse_nonsym examples/reverse/nonsym/*.cc)
examples(reverse_nonsym)
file(GLOB reverse_sym examples/reverse/sym/*.cc)
examples(reverse_sym)

# examples band
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/band/)
file(GLOB band_complex examples/band/complex/*.cc)
examples(band_complex)
file(GLOB band_nonsym examples/band/nonsym/*.cc)
examples(band_nonsym)
file(GLOB band_sym examples/band/sym/*.cc)
examples(band_sym)

# examples dense
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/dense/)
file(GLOB dense_complex examples/dense/complex/*.cc)
examples(dense_complex)
file(GLOB dense_nonsym examples/dense/nonsym/*.cc)
examples(dense_nonsym)
file(GLOB dense_sym examples/dense/sym/*.cc)
examples(dense_sym)




if (SUPERLU)

  # examples areig
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/areig/)
  file(GLOB areig_complex examples/areig/complex/*.cc)
  examplesslu(areig_complex)
  file(GLOB areig_nonsym examples/areig/nonsym/*.cc)
  examplesslu(areig_nonsym)
  file(GLOB areig_sym examples/areig/sym/*.cc)
  examplesslu(areig_sym)

  # examples superlu
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/superlu/)
  file(GLOB superlu_complex examples/superlu/complex/*.cc)
  examplesslu(superlu_complex)
  file(GLOB superlu_nonsym examples/superlu/nonsym/*.cc)
  examplesslu(superlu_nonsym)
  file(GLOB superlu_sym examples/superlu/sym/*.cc)
  examplesslu(superlu_sym)

  # examples harwell (needs SuperLU)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/harwell/)
  file(GLOB harwell_complex examples/harwell/complex/*.cc)
  examplesslu(harwell_complex)
  file(GLOB harwell_nonsym examples/harwell/nonsym/*.cc)
  examplesslu(harwell_nonsym)
  file(GLOB harwell_sym examples/harwell/sym/*.cc)
  examplesslu(harwell_sym)
  # also copy binaries
  add_custom_target(harwellbin
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/examples/harwell/complex/mhd1280a.cua.gz" ${CMAKE_BINARY_DIR}/examples/harwell/ 
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/examples/harwell/complex/mhd1280b.cua.gz" ${CMAKE_BINARY_DIR}/examples/harwell/ 
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/examples/harwell/nonsym/mhd416a.rua" ${CMAKE_BINARY_DIR}/examples/harwell/ 
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/examples/harwell/nonsym/mhd416b.rua" ${CMAKE_BINARY_DIR}/examples/harwell/ 
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/examples/harwell/sym/lund_a.rsa" ${CMAKE_BINARY_DIR}/examples/harwell/ 
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/examples/harwell/sym/lund_b.rsa" ${CMAKE_BINARY_DIR}/examples/harwell/ 
    COMMENT "Copying Harwell binary examples"
    SOURCES examples/harwell/complex/mhd1280a.cua.gz 
      examples/harwell/complex/mhd1280b.cua.gz 
      examples/harwell/nonsym/mhd416a.rua
      examples/harwell/nonsym/mhd416b.rua
      examples/harwell/sym/lund_a.rsa
      examples/harwell/sym/lund_b.rsa
  )
  add_dependencies(examples harwellbin)
endif()

if (UMFPACK)
  
  # examples umfpack
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/umfpack/)
  file(GLOB umfpack_complex examples/umfpack/complex/*.cc)
  #examplesls(umfpack_complex)
  file(GLOB umfpack_nonsym examples/umfpack/nonsym/*.cc)
  #examplesls(umfpack_nonsym)
  file(GLOB umfpack_sym examples/umfpack/sym/*.cc)
  examplesls(umfpack_sym)
    
endif()

if (CHOLMOD)

  # examples cholmod
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples/cholmod/)
  file(GLOB cholmod_sym examples/cholmod/sym/*.cc)
  examplesls(cholmod_sym)

endif()

add_custom_target( Message ALL
    COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "Run either \\\"make examples\\\" or \\\"make install\\\"."
    COMMENT "Header Library: Nothing to compile!" )

add_library(arpackpp INTERFACE)
target_include_directories(arpackpp INTERFACE include)

add_custom_target (examples)

set(COMMON_INCLUDES
  ${CMAKE_CURRENT_SOURCE_DIR}/matrices/complex
  ${CMAKE_CURRENT_SOURCE_DIR}/matrices/nonsym
  ${CMAKE_CURRENT_SOURCE_DIR}/matrices/sym
  ${CMAKE_CURRENT_SOURCE_DIR}/matprod
  ${CMAKE_CURRENT_SOURCE_DIR}/matprod/complex
  ${CMAKE_CURRENT_SOURCE_DIR}/matprod/nonsym
  ${CMAKE_CURRENT_SOURCE_DIR}/matprod/sym
  ${CMAKE_CURRENT_SOURCE_DIR}/areig
)

function(examples list_name)
  foreach(l ${${list_name}})
    get_filename_component(lwe ${l} NAME_WE)
    add_executable(${lwe} ${l})
    target_include_directories(${lwe} PRIVATE ${COMMON_INCLUDES})
    target_link_libraries(${lwe}
      PUBLIC
        ${ARPACK_LIB}
        arpackpp
        $<BUILD_INTERFACE:LAPACK::LAPACK>
        $<BUILD_INTERFACE:BLAS::BLAS>)
    add_dependencies (examples ${lwe})
  endforeach()
endfunction(examples)

function(examplesslu list_name)
  foreach(l ${${list_name}})
    get_filename_component(lwe ${l} NAME_WE)
    add_executable(${lwe} ${l})
    target_include_directories(${lwe} PRIVATE ${COMMON_INCLUDES})
    target_link_libraries(${lwe}
      PUBLIC
        ${ARPACK_LIB}
        arpackpp
        ${SUPERLU_LIB}
        $<BUILD_INTERFACE:LAPACK::LAPACK>
        $<BUILD_INTERFACE:BLAS::BLAS>)
    add_dependencies (examples ${lwe})
  endforeach()
endfunction(examplesslu)

function(examplesls list_name)
  foreach(l ${${list_name}})
    get_filename_component(lwe ${l} NAME_WE)
    add_executable(${lwe} ${l})
    target_include_directories(${lwe} PRIVATE ${COMMON_INCLUDES})
    target_link_libraries(${lwe}
      PUBLIC
        ${ARPACK_LIB}
        arpackpp
        ${SSPARSE_LIBS}
        $<BUILD_INTERFACE:LAPACK::LAPACK>
        $<BUILD_INTERFACE:BLAS::BLAS>)
    add_dependencies (examples ${lwe})
  endforeach()
endfunction(examplesls)

# examples product
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/product/)
file(GLOB product_complex product/complex/*.cc)
examples(product_complex)
file(GLOB product_nonsym product/nonsym/*.cc)
examples(product_nonsym)
file(GLOB product_simple product/simple/*.cc)
examples(product_simple)
file(GLOB product_sym product/sym/*.cc)
examples(product_sym)

# examples reverse
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/reverse/)
file(GLOB reverse_complex reverse/complex/*.cc)
examples(reverse_complex)
file(GLOB reverse_nonsym reverse/nonsym/*.cc)
examples(reverse_nonsym)
file(GLOB reverse_sym reverse/sym/*.cc)
examples(reverse_sym)

# examples band
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/band/)
file(GLOB band_complex band/complex/*.cc)
examples(band_complex)
file(GLOB band_nonsym band/nonsym/*.cc)
examples(band_nonsym)
file(GLOB band_sym band/sym/*.cc)
examples(band_sym)

# examples dense
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/dense/)
file(GLOB dense_complex dense/complex/*.cc)
examples(dense_complex)
file(GLOB dense_nonsym dense/nonsym/*.cc)
examples(dense_nonsym)
file(GLOB dense_sym dense/sym/*.cc)
examples(dense_sym)




if (ENABLE_SUPERLU)

  # examples areig
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/areig/)
  file(GLOB areig_complex areig/complex/*.cc)
  examplesslu(areig_complex)
  file(GLOB areig_nonsym areig/nonsym/*.cc)
  examplesslu(areig_nonsym)
  file(GLOB areig_sym areig/sym/*.cc)
  examplesslu(areig_sym)

  # examples superlu
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/superlu/)
  file(GLOB superlu_complex superlu/complex/*.cc)
  examplesslu(superlu_complex)
  file(GLOB superlu_nonsym superlu/nonsym/*.cc)
  examplesslu(superlu_nonsym)
  file(GLOB superlu_sym superlu/sym/*.cc)
  examplesslu(superlu_sym)

  # examples harwell (needs SuperLU)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/harwell/)
  file(GLOB harwell_complex harwell/complex/*.cc)
  examplesslu(harwell_complex)
  file(GLOB harwell_nonsym harwell/nonsym/*.cc)
  examplesslu(harwell_nonsym)
  file(GLOB harwell_sym harwell/sym/*.cc)
  examplesslu(harwell_sym)
  # also copy binaries
  add_custom_target(harwellbin
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/harwell/complex/mhd1280a.cua.gz" ${CMAKE_BINARY_DIR}/harwell/ 
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/harwell/complex/mhd1280b.cua.gz" ${CMAKE_BINARY_DIR}/harwell/ 
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/harwell/nonsym/mhd416a.rua" ${CMAKE_BINARY_DIR}/harwell/ 
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/harwell/nonsym/mhd416b.rua" ${CMAKE_BINARY_DIR}/harwell/ 
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/harwell/sym/lund_a.rsa" ${CMAKE_BINARY_DIR}/harwell/ 
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/harwell/sym/lund_b.rsa" ${CMAKE_BINARY_DIR}/harwell/ 
    COMMENT "Copying Harwell binary examples"
    SOURCES harwell/complex/mhd1280a.cua.gz 
      harwell/complex/mhd1280b.cua.gz 
      harwell/nonsym/mhd416a.rua
      harwell/nonsym/mhd416b.rua
      harwell/sym/lund_a.rsa
      harwell/sym/lund_b.rsa
  )
  add_dependencies(examples harwellbin)
endif()

if (ENABLE_UMFPACK)
  
  # examples umfpack
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/umfpack/)
  file(GLOB umfpack_complex umfpack/complex/*.cc)
  #examplesls(umfpack_complex)
  file(GLOB umfpack_nonsym umfpack/nonsym/*.cc)
  #examplesls(umfpack_nonsym)
  file(GLOB umfpack_sym umfpack/sym/*.cc)
  examplesls(umfpack_sym)
    
endif()

if (ENABLE_CHOLMOD)

  # examples cholmod
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/cholmod/)
  file(GLOB cholmod_sym cholmod/sym/*.cc)
  examplesls(cholmod_sym)

endif()

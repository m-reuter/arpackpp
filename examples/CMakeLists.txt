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
  # Get optional solver argument
  set (solver ${ARGN})

  foreach(l ${${list_name}})
    get_filename_component(lwe ${l} NAME_WE)
    add_executable(${lwe} ${l})
    target_link_libraries(${lwe}
      PRIVATE
        $<BUILD_INTERFACE:LAPACK::LAPACK>
        $<BUILD_INTERFACE:BLAS::BLAS>
        arpackpp)

    if (solver STREQUAL "superlu")
      target_link_libraries(${lwe} PRIVATE $<BUILD_INTERFACE:superlu::superlu>)
    elseif (solver STREQUAL "cholmod")
      target_link_libraries(${lwe} PRIVATE $<BUILD_INTERFACE:$<IF:${suitesparse_static},SuiteSparse::CHOLMOD_static,SuiteSparse::CHOLMOD>>)
      target_include_directories(${lwe} PRIVATE ${CHOLMOD_INCLUDE_DIR})
    elseif (solver STREQUAL "umfpack")
      target_link_libraries(${lwe} PRIVATE $<BUILD_INTERFACE:$<IF:${suitesparse_static},SuiteSparse::UMFPACK_static,SuiteSparse::UMFPACK>>)
    endif ()

    target_include_directories(${lwe} PRIVATE ${COMMON_INCLUDES})
    add_dependencies (examples ${lwe})
  endforeach()
endfunction()

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
  examples(areig_complex "superlu")
  file(GLOB areig_nonsym areig/nonsym/*.cc)
  examples(areig_nonsym "superlu")
  file(GLOB areig_sym areig/sym/*.cc)
  examples(areig_sym "superlu")

  # examples superlu
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/superlu/)
  file(GLOB superlu_complex superlu/complex/*.cc)
  examples(superlu_complex "superlu")
  file(GLOB superlu_nonsym superlu/nonsym/*.cc)
  examples(superlu_nonsym "superlu")
  file(GLOB superlu_sym superlu/sym/*.cc)
  examples(superlu_sym "superlu")

  # examples harwell (needs SuperLU)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/harwell/)
  file(GLOB harwell_complex harwell/complex/*.cc)
  examples(harwell_complex "superlu")
  file(GLOB harwell_nonsym harwell/nonsym/*.cc)
  examples(harwell_nonsym "superlu")
  file(GLOB harwell_sym harwell/sym/*.cc)
  examples(harwell_sym "superlu")
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
  #examples(umfpack_complex "umfpack")
  file(GLOB umfpack_nonsym umfpack/nonsym/*.cc)
  #examples(umfpack_nonsym "umfpack")
  file(GLOB umfpack_sym umfpack/sym/*.cc)
  examples(umfpack_sym "umfpack")
    
endif()

if (ENABLE_CHOLMOD)

  # examples cholmod
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/cholmod/)
  file(GLOB cholmod_sym cholmod/sym/*.cc)
  examples(cholmod_sym "cholmod")

endif()

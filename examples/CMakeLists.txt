add_custom_target (tests)

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

function(setup_tests list)
  # Optional arguments:
  # - SOLVER = name of solver (string)
  # - INPUT  = arguments passed to test executable (list of strings)
  cmake_parse_arguments(test "" "SOLVER" "INPUT" ${ARGN} )

  foreach(file ${${list}})
    get_filename_component(target ${file} NAME_WE)
    add_executable(${target} ${file})
    target_link_libraries(${target}
      PRIVATE
        $<BUILD_INTERFACE:LAPACK::LAPACK>
        $<BUILD_INTERFACE:BLAS::BLAS>
        arpackpp)

    if (test_SOLVER STREQUAL "superlu")
      target_link_libraries(${target} PRIVATE $<BUILD_INTERFACE:superlu::superlu>)
    elseif (test_SOLVER STREQUAL "cholmod")
      target_link_libraries(${target} PRIVATE $<BUILD_INTERFACE:$<IF:$<BOOL:${ENABLE_SUITESPARSE_STATIC}>,SuiteSparse::CHOLMOD_static,SuiteSparse::CHOLMOD>>)
      target_include_directories(${target} PRIVATE ${CHOLMOD_INCLUDE_DIR})
    elseif (test_SOLVER STREQUAL "umfpack")
      target_link_libraries(${target} PRIVATE $<BUILD_INTERFACE:$<IF:$<BOOL:${ENABLE_SUITESPARSE_STATIC}>,SuiteSparse::UMFPACK_static,SuiteSparse::UMFPACK>>)
      target_include_directories(${target} PRIVATE ${UMFPACK_INCLUDE_DIR})
    endif ()

    target_include_directories(${target} PRIVATE ${COMMON_INCLUDES})
    add_dependencies (tests ${target})

    add_test(NAME ${target}_test
      COMMAND ${target} ${test_INPUT}
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/data")

  endforeach()
endfunction()

# examples product
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/product/)
file(GLOB product_complex product/complex/*.cc)
setup_tests(product_complex)
file(GLOB product_nonsym product/nonsym/*.cc)
setup_tests(product_nonsym)
file(GLOB product_simple product/simple/*.cc)
setup_tests(product_simple)
file(GLOB product_sym product/sym/*.cc)
setup_tests(product_sym)

# examples reverse
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/reverse/)
file(GLOB reverse_complex reverse/complex/*.cc)
setup_tests(reverse_complex)
file(GLOB reverse_nonsym reverse/nonsym/*.cc)
setup_tests(reverse_nonsym)
file(GLOB reverse_sym reverse/sym/*.cc)
setup_tests(reverse_sym)

# examples band
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/band/)
file(GLOB band_complex band/complex/*.cc)
setup_tests(band_complex)
file(GLOB band_nonsym band/nonsym/*.cc)
setup_tests(band_nonsym)
file(GLOB band_sym band/sym/*.cc)
setup_tests(band_sym)

# examples dense
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/dense/)
file(GLOB dense_complex dense/complex/*.cc)
setup_tests(dense_complex)
file(GLOB dense_nonsym dense/nonsym/*.cc)
setup_tests(dense_nonsym)
file(GLOB dense_sym dense/sym/*.cc)
setup_tests(dense_sym)

# copy test data
file(ARCHIVE_EXTRACT INPUT "${CMAKE_CURRENT_SOURCE_DIR}/dense/nonsym/matrix.zip" DESTINATION "${CMAKE_BINARY_DIR}/data")

if (ENABLE_SUPERLU)

  # examples areig
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/areig/)
  file(GLOB areig_complex areig/complex/*.cc)
  setup_tests(areig_complex SOLVER "superlu")
  file(GLOB areig_nonsym areig/nonsym/*.cc)
  setup_tests(areig_nonsym SOLVER "superlu")
  file(GLOB areig_sym areig/sym/*.cc)
  setup_tests(areig_sym SOLVER "superlu")

  # examples superlu
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/superlu/)
  file(GLOB superlu_complex superlu/complex/*.cc)
  setup_tests(superlu_complex SOLVER "superlu")
  file(GLOB superlu_nonsym superlu/nonsym/*.cc)
  setup_tests(superlu_nonsym SOLVER "superlu")
  file(GLOB superlu_sym superlu/sym/*.cc)
  setup_tests(superlu_sym SOLVER "superlu")

  # examples harwell (needs SuperLU)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/harwell/)

  set(file "harwell/complex/hcompstd.cc")
  set(args "-n" "4" "mhd1280a.cua")
  setup_tests(file SOLVER "superlu" INPUT ${args})
  
  set(file "harwell/complex/hcompgen.cc")
  set(args "-n" "4" "mhd1280a.cua" "mhd1280b.cua")
  setup_tests(file SOLVER "superlu" INPUT ${args})

  set(file "harwell/nonsym/hnsymstd.cc")
  set(args "-n" "4" "mhd416a.rua")
  setup_tests(file SOLVER "superlu" INPUT ${args})
  
  set(file "harwell/nonsym/hnsymgen.cc")
  set(args "-n" "4" "mhd416a.rua" "mhd416b.rua")
  setup_tests(file SOLVER "superlu" INPUT ${args})
  
  set(file "harwell/sym/hsymstd.cc")
  set(args "-n" "4" "lund_a.rsa")
  setup_tests(file SOLVER "superlu" INPUT ${args})
  
  set(file "harwell/sym/hsymgen.cc")
  set(args "-n" "4" "lund_a.rsa" "lund_b.rsa")
  setup_tests(file SOLVER "superlu" INPUT ${args})

  # copy test data
  file(ARCHIVE_EXTRACT INPUT "${CMAKE_CURRENT_SOURCE_DIR}/harwell/complex/mhd1280a.zip" DESTINATION "${CMAKE_BINARY_DIR}/data")
  file(ARCHIVE_EXTRACT INPUT "${CMAKE_CURRENT_SOURCE_DIR}/harwell/complex/mhd1280b.zip" DESTINATION "${CMAKE_BINARY_DIR}/data")
  file(COPY
    "${CMAKE_CURRENT_SOURCE_DIR}/harwell/nonsym/mhd416a.rua"
    "${CMAKE_CURRENT_SOURCE_DIR}/harwell/nonsym/mhd416b.rua"
    "${CMAKE_CURRENT_SOURCE_DIR}/harwell/sym/lund_a.rsa"
    "${CMAKE_CURRENT_SOURCE_DIR}/harwell/sym/lund_b.rsa"
  DESTINATION "${CMAKE_BINARY_DIR}/data")
endif()

if (ENABLE_UMFPACK)
  
  # examples umfpack
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/umfpack/)
  file(GLOB umfpack_complex umfpack/complex/*.cc)
  setup_tests(umfpack_complex SOLVER "umfpack")
  file(GLOB umfpack_nonsym umfpack/nonsym/*.cc)
  setup_tests(umfpack_nonsym SOLVER "umfpack")
  file(GLOB umfpack_sym umfpack/sym/*.cc)
  setup_tests(umfpack_sym SOLVER "umfpack")
    
endif()

if (ENABLE_CHOLMOD)

  # examples cholmod
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/cholmod/)
  file(GLOB cholmod_sym cholmod/sym/*.cc)
  setup_tests(cholmod_sym SOLVER "cholmod")

endif()

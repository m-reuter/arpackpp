cmake_minimum_required (VERSION 3.18)

project(arpackpp VERSION 2.4.0
                 DESCRIPTION "ARPACK++: C++ interface to ARPACK"
                 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(arpackpp_MAJOR_VERSION 2)
set(arpackpp_MINOR_VERSION 4)
set(arpackpp_PATCH_VERSION 0)
set(arpackpp_VERSION ${arpackpp_MAJOR_VERSION}.${arpackpp_MINOR_VERSION}.${arpackpp_PATCH_VERSION})

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/external/lib/libumfpack.a" AND
   EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/external/lib/libcholmod.a")
  set (suitesparse_static ON)
else()
  set (suitesparse_static OFF)
endif()

option(ENABLE_TESTS "Build tests (examples)" ON)
option(ENABLE_FORTRAN "Enable Fortran language (for static linking of ARPACK)" OFF)
option(ENABLE_SUPERLU "Enable SUPERLU" OFF)
option(ENABLE_UMFPACK "Enable UMFPACK" OFF)
option(ENABLE_CHOLMOD "Enable CHOLMOD" OFF)
option(ENABLE_SUITESPARSE_STATIC "Enable linking SuiteSparse static targets" ${suitesparse_static})
option(INSTALL_ARPACKPP_CMAKE_TARGET "Enable the creation of CMake config targets" ON)

if (ENABLE_FORTRAN)
  enable_language(Fortran)
endif()

# In case the dependencies are installed locally using the shell scripts, we need to
# let CMake know, where to look for the config-file packages.
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_SOURCE_DIR}/external)

# In case dependencies were installed by a package manager, we need to let CMake know
# where to look for the modules.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

# Custom macro to find packages. First try config mode. If the target isn't found,
# use module mode.
macro(find_package_custom name target)
  find_package(${name} CONFIG QUIET)
  
  if (NOT TARGET ${target})
    find_package(${name} REQUIRED)
  else()
    get_target_property(location ${target} LOCATION)
    message( STATUS "Found ${name} [config]: ${location}" )
  endif()
endmacro()

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
  if (ENABLE_SUITESPARSE_STATIC)
    message(STATUS "Linking SuiteSparse static")
    find_package_custom ( SuiteSparse_config SuiteSparse::SuiteSparseConfig_static )
    find_package_custom ( AMD SuiteSparse::AMD_static )
    find_package_custom ( CAMD SuiteSparse::CAMD_static )
    find_package_custom ( CCOLAMD SuiteSparse::CCOLAMD_static )
    find_package_custom ( COLAMD SuiteSparse::COLAMD_static )
    find_package_custom ( CHOLMOD SuiteSparse::CHOLMOD_static )

    if (NOT TARGET SuiteSparse::CHOLMOD_static)
      message(FATAL_ERROR "CHOLMOD_static not found")
    endif()
  else()
    find_package_custom ( SuiteSparse_config SuiteSparse::SuiteSparseConfig )
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
  if (ENABLE_SUITESPARSE_STATIC)
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

include(GNUInstallDirs)

add_library(arpackpp INTERFACE)
add_library(arpackpp::arpackpp ALIAS arpackpp)

target_link_libraries(arpackpp INTERFACE ARPACK::ARPACK)

# Adding the install interface generator expression makes sure that the include
# files are installed to the proper location (provided by GNUInstallDirs)
target_include_directories(arpackpp INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>/arpackpp)

target_compile_features(arpackpp INTERFACE cxx_std_17)

# Examples

if(ENABLE_TESTS)
  enable_testing()
  add_subdirectory(examples)
endif()

# Install

install (TARGETS arpackpp EXPORT arpackppTargets)

if(INSTALL_ARPACKPP_CMAKE_TARGET)
  include(CMakePackageConfigHelpers)

  set(ARPACKPP_CONFIG_INSTALL_DIR ${CMAKE_INSTALL_DATADIR}/arpackpp/cmake)

  configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/arpackppConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/arpackppConfig.cmake
    INSTALL_DESTINATION ${ARPACKPP_CONFIG_INSTALL_DIR}
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
  )

  write_basic_package_version_file(arpackppConfigVersion.cmake
    VERSION ${arpackpp_VERSION}
    COMPATIBILITY SameMajorVersion
    ARCH_INDEPENDENT)

  export(TARGETS arpackpp NAMESPACE arpackpp:: FILE arpackppTargets.cmake)

  install(EXPORT arpackppTargets NAMESPACE arpackpp:: DESTINATION ${ARPACKPP_CONFIG_INSTALL_DIR})
  install (FILES
    ${CMAKE_CURRENT_BINARY_DIR}/arpackppConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/arpackppConfigVersion.cmake
    DESTINATION ${ARPACKPP_CONFIG_INSTALL_DIR})
endif()

file (GLOB HEADERS "include/*.h")
install (FILES ${HEADERS} DESTINATION include/arpackpp)

file (GLOB_RECURSE EXAMPLE_HEADERS "examples/*.h")
install (FILES ${EXAMPLE_HEADERS} DESTINATION include/arpackpp)

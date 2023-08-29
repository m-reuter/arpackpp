find_package(PkgConfig QUIET)

if (PKG_CONFIG_FOUND)
  pkg_check_modules(PKG_LIBARPACK QUIET arpack)
endif()

find_path(ARPACK_INCLUDE_DIRS arpack.h
  ${PKG_LIBARPACK_INCLUDE_DIRS}
  /usr/include
  /usr/local/include
)

find_library(ARPACK_LIBRARIES
  NAMES
    arpack libarpack
  PATHS
    ${PKG_LIBARPACK_LIBRARY_DIRS}
    /usr/lib
    /usr/local/lib
)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(ARPACK
    FOUND_VAR
        ARPACK_FOUND
    REQUIRED_VARS
        ARPACK_LIBRARIES
    #    ARPACK_INCLUDE_DIRS
    #VERSION_VAR
    #    ARPACK_VERSION
)

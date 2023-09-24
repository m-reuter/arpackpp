if (NOT WIN32)
  find_package(PkgConfig QUIET)
  if (PKG_CONFIG_FOUND)
    pkg_check_modules(SUPERLU superlu QUIET)
  endif ()
endif (NOT WIN32)

if ( NOT SUPERLU_FOUND )
  find_path(SUPERLU_INCLUDE_DIRS
    NAMES slu_util.h
    PATHS
      /usr/include
      /usr/local/include
    PATH_SUFFIXES
      superlu)

  find_library(SUPERLU_LIBRARIES
    NAMES superlu
    PATHS
      /usr/lib
      /usr/local/lib)

  if ( NOT "${SUPERLU_LIBRARIES}" STREQUAL "")
    set (SUPERLU_FOUND ${SUPERLU_LIBRARIES_FOUND})
    set (SUPERLU_LINK_LIBRARIES ${SUPERLU_LINK_LIBRARIES} ${SUPERLU_LIBRARIES})
  endif ()
endif ()

include (FindPackageHandleStandardArgs)

find_package_handle_standard_args ( superlu
    REQUIRED_VARS SUPERLU_LINK_LIBRARIES SUPERLU_INCLUDE_DIRS
    VERSION_VAR SUPERLU_VERSION
)

mark_as_advanced (
    SUPERLU_INCLUDE_DIRS
    SUPERLU_LIBRARIES
)

if ( SUPERLU_FOUND AND NOT TARGET superlu::superlu )
  add_library(superlu::superlu INTERFACE IMPORTED GLOBAL)
  set_target_properties(superlu::superlu
    PROPERTIES
      LOCATION "${SUPERLU_LIBRARIES}"
      INTERFACE_INCLUDE_DIRECTORIES "${SUPERLU_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "${SUPERLU_LINK_LIBRARIES}")
endif ()


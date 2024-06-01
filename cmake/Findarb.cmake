
set(ARB_ROOT_DIR "${ARB_ROOT_DIR}"  CACHE PATH "Directory to search for arb" )

find_package(PkgConfig QUIET)
if( PkgConfig_FOUND )
  pkg_search_module(PC_ARB QUIET arb flint)
endif()

find_path( ARB_INCLUDE_DIR
  NAMES arb.h
  PATHS "${ARB_ROOT_DIR}"
  HINTS ${PC_ARB_INCLUDEDIR} ${PC_ARB_INCLUDE_DIRS}
  PATH_SUFFIXES flint arb
  )
find_library( ARB_LIBRARY
  NAMES arb flint
  PATHS "${ARB_ROOT_DIR}"
  HINTS ${PC_ARB_LIBDIR} ${PC_ARB_LIBRARY_DIRS}
  )

set ( _VERSION_FILE ${ARB_INCLUDE_DIR}/arb.h )
if ( EXISTS ${_VERSION_FILE} )
  file ( STRINGS ${_VERSION_FILE} _VERSION_LINE REGEX "define[ ]+ARB_VERSION" )
  if ( _VERSION_LINE )
    string ( REGEX REPLACE ".*define[ ]+ARB_VERSION[ ]+\"([^\"]*)\".*" "\\1" ARB_VERSION "${_VERSION_LINE}" )
  endif ()
endif ()
unset ( _VERSION_FILE )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args( arb
  FOUND_VAR ARB_FOUND
  REQUIRED_VARS
    ARB_LIBRARY
    ARB_INCLUDE_DIR
  VERSION_VAR ARB_VERSION
  )

if(ARB_FOUND)
  set(ARB_INCLUDE_DIRS ${ARB_INCLUDE_DIR})
  set(ARB_LIBRARIES ${ARB_LIBRARY})
  if(NOT TARGET arb::arb)
    add_library(arb::arb UNKNOWN IMPORTED)
    set_target_properties( arb::arb
      PROPERTIES
        IMPORTED_LOCATION "${ARB_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${ARB_INCLUDE_DIR}"
      )
  endif()
  mark_as_advanced(ARB_ROOT_DIR)
endif()

mark_as_advanced(ARB_INCLUDE_DIR ARB_LIBRARY)

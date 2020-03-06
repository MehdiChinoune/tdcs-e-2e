
option( Profiling "Enable profiling with gprof" OFF )
if( Profiling )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -p" )
endif()

set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules )

file( GLOB srcs CONFIGURE_DEPENDS *.F90 )
list( FILTER srcs EXCLUDE REGEX "main.F90"  )

add_library( tdcs_e_2e SHARED ${srcs} )

option( OpenMP "Enable OpenMP" OFF )
if( OpenMP )
  find_package( OpenMP REQUIRED )
  target_link_libraries( tdcs_e_2e PUBLIC OpenMP::OpenMP_Fortran )
endif()

add_executable( main "main.F90" )
target_link_libraries( main tdcs_e_2e )

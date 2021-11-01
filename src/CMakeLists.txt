
file( GLOB srcs CONFIGURE_DEPENDS *.F90 )
list( FILTER srcs EXCLUDE REGEX "main.F90"  )

add_library( tdcs_e_2e SHARED ${srcs} )

if( OpenMP )
  target_link_libraries( tdcs_e_2e PUBLIC OpenMP::OpenMP_Fortran )
endif()

add_executable( main "main.F90" )
target_link_libraries( main tdcs_e_2e )

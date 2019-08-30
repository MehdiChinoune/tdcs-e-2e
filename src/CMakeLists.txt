

if( OpenMP )
  find_package( OpenMP )
  if( OpenMP_Fortran_FOUND )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}" )
  elseif( MINGW )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp" )
  endif()
endif()

set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules )

file( GLOB srcs CONFIGURE_DEPENDS *.F90 )
list( FILTER srcs EXCLUDE REGEX "main.F90"  )

add_library( tdcs_e_2e SHARED ${srcs} )

add_executable( main "main.F90" )
target_link_libraries( main tdcs_e_2e )


set(tests test_conhyp_01 )

INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR}/modules )

foreach( tst ${tests} )
  add_executable( ${tst} "${tst}.F90" )
  target_link_libraries( ${tst} tdcs_e_2e )
  if( ENABLE_OPENMP )
    target_link_libraries( ${tst} OpenMP::OpenMP_Fortran )
  endif()
  add_test( NAME ${tst} COMMAND ./${tst} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  set_tests_properties( ${tst} PROPERTIES FAIL_REGULAR_EXPRESSION "[^a-z]Failed;FAILED" )
endforeach()

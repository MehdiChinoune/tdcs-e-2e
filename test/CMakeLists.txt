
find_package(arb 2.20 REQUIRED)

add_library(arb_bindings OBJECT arb_functions.F90)
target_link_libraries(arb_bindings arb::arb tdcs_e_2e)

set(tests
  test_conhyp_01
  test_conhyp_opt
  test_arb_hypgeom_01
  test_arb_hypgeom_02
  test_arb_hypgeom_03
)

INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR}/modules )

foreach( tst ${tests} )
  add_executable( ${tst} "${tst}.F90" )
  target_link_libraries( ${tst} tdcs_e_2e arb_bindings )
  if( ENABLE_OPENMP )
    target_link_libraries( ${tst} OpenMP::OpenMP_Fortran )
  endif()
  add_test( NAME ${tst} COMMAND ./${tst} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  set_tests_properties( ${tst} PROPERTIES FAIL_REGULAR_EXPRESSION "[^a-z]Failed;FAILED" )
endforeach()

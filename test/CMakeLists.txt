set(tests
  test_adaptive_integration
  test_orbit_wf
)

if(BUILD_SLOW_TESTS)
  list(APPEND tests test_eikr_element)
endif()

if(NOT (WIN32 AND "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "LLVMFlang"))
  find_package(arb 2.20 REQUIRED)

  add_library(arb_bindings OBJECT arb_functions.F90)
  target_link_libraries(arb_bindings arb::arb tdcs_e_2e)

  list(APPEND tests
    test_conhyp_01
    test_conhyp_opt
  )
  if(BUILD_SLOW_TESTS)
    list(APPEND tests
      test_arb_hypgeom_01
      test_arb_hypgeom_02
      test_arb_hypgeom_03
    )
  endif()
endif()

INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR}/modules )

foreach( tst ${tests} )
  add_executable( ${tst} "${tst}.F90" )
  target_link_libraries( ${tst} tdcs_e_2e )
  if(NOT (WIN32 AND "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "LLVMFlang"))
    target_link_libraries( ${tst} arb_bindings )
  endif()
  if( ENABLE_OPENMP )
    target_link_libraries( ${tst} OpenMP::OpenMP_Fortran )
  endif()
  add_test( NAME ${tst} COMMAND ./${tst} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  set_tests_properties( ${tst} PROPERTIES FAIL_REGULAR_EXPRESSION "[^a-z]Failed;FAILED" )
endforeach()

if(TARGET test_eikr_element AND WIN32)
  if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "LLVMFlang")
    target_link_options(test_eikr_element PRIVATE "-Wl,--stack,134217728")
  elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
    target_link_options(test_eikr_element PRIVATE "/STACK:536870912")
  endif()
endif()

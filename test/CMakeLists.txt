
FILE(GLOB tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS test*.F90 )

INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR}/modules )

foreach( tst ${tests} )
  STRING(REPLACE ".F90" "" exe ${tst})
  add_executable( ${exe} ${tst} )
  target_link_libraries( ${exe} tdcs_e_2e )
  add_test( NAME ${exe} COMMAND ./${exe} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  set_tests_properties( ${exe} PROPERTIES FAIL_REGULAR_EXPRESSION "[^a-z]Failed;FAILED" )
endforeach()

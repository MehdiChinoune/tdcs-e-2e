
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}" )

set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/obj/${CMAKE_BUILD_TYPE})

set(srcs main.F90 constants.F90 trigo.F90 special_functions.F90 integration.F90 input.F90 utils.F90 fdcs_e2e.F90 )

add_executable(main ${srcs})

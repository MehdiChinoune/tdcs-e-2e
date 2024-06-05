
set( srcs
  constants.F90
  conhyp_h.F90
  conhyp.F90
  input_h.F90
  input.F90
  integration_h.F90
  integration.F90
  trigo_h.F90
  trigo.F90
  special_functions_h.F90
  special_functions.F90
  utils_h.F90
  utils.F90
  fdcs_e2e_h.F90
  fdcs_e2e.F90
  tcs_e2e_h.F90
  tcs_e2e.F90
  types.F90
)

if(WIN32 AND "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "LLVMFlang")
  add_library( tdcs_e_2e STATIC ${srcs} )
else()
  add_library( tdcs_e_2e SHARED ${srcs} )
endif()

if( TDCS_OPENMP )
  target_link_libraries( tdcs_e_2e PUBLIC OpenMP::OpenMP_Fortran )
endif()

add_executable( main "main.F90" )
target_link_libraries( main tdcs_e_2e )

if(WIN32)
  if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "LLVMFlang")
    target_link_options(main PRIVATE "-Wl,--stack,134217728")
  elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
    target_link_options(main PRIVATE "/STACK:536870912")
  endif()
endif()

cmake_minimum_required(VERSION 3.7)

project(TDCS-e-2e Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set( CMAKE_BUILD_TYPE "Release" )
  message(STATUS "No CMAKE_BUILD_TYPE selected, default to ${CMAKE_BUILD_TYPE}")
endif()

if ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" )

  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wimplicit-procedure -Wunused-parameter -Wunused-but-set-parameter -Wunused-but-set-variable -std=f2018" )
  set( CMAKE_Fortran_FLAGS_FAST "-Ofast -flto" )
  set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -flto" )
  set( CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fcheck=all -ffpe-trap=zero,overflow,invalid -fbacktrace" )

elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel" )

  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -warn all -stand f18 -standard-semantics" )
  set( CMAKE_Fortran_FLAGS_FAST "-Ofast -ipo" )
  set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -ipo" )
  set( CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -check all -fpe-all=0 -traceback" )

elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "Flang" )

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -Wall -Mstandard -O3" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -Wall -Mstandard -g" )

elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "PGI" )

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -Mstandard -Mallocatable=03 -fast" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -Mstandard -Mallocatable=03 -g" )

endif ()

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}" )

add_subdirectory(src)

file(COPY Data DESTINATION ./ )
file(COPY input_example.dat DESTINATION ./ )

enable_testing()
add_subdirectory(test)

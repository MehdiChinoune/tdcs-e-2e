cmake_minimum_required(VERSION 3.7)

project(TDCS-e-2e Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set( CMAKE_BUILD_TYPE "Release" )
  message(STATUS "No CMAKE_BUILD_TYPE selected, default to ${CMAKE_BUILD_TYPE}")
endif()

Option( TARGET "Generate target optimized binaries" ON )

if ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" )

  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wimplicit-procedure -Wunused-parameter -std=f2018" )
  set( CMAKE_Fortran_FLAGS_FAST "-Ofast -flto" )
  set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -flto" )
  set( CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fcheck=all -ffpe-trap=zero,overflow,invalid -fbacktrace" )
  if( TARGET )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native" )
  endif()

elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel" )

  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -warn all -stand f18 -standard-semantics" )
  set( CMAKE_Fortran_FLAGS_FAST "-Ofast -ipo" )
  set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -ipo" )
  set( CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -check all -fpe-all=0 -traceback" )
  if( TARGET )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xHost" )
  endif()

else ()

  set( CMAKE_Fortran_FLAGS_RELEASE "-O3" )
  set( CMAKE_Fortran_FLAGS_DEBUG "-O0 -g" )

endif ()

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}" )

add_subdirectory(src)

file(COPY data DESTINATION ./ )
file(COPY input_example.dat DESTINATION ./ )

enable_testing()
add_subdirectory(test)

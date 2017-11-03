cmake_minimum_required(VERSION 3.7)

project(TDCS-e-2e Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set( CMAKE_BUILD_TYPE "Release" )
  message(STATUS "No CMAKE_BUILD_TYPE selected, default to ${CMAKE_BUILD_TYPE}")
endif()

if ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" )

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -Wall -O3" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -Wall -g -fcheck=bounds" )
  
elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "INTEL" )

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -warn all -Ofast -ipo -xHost" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -warn all -g" )
  
elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "FLANG" )

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -Wall -Mstandard -Mallocatable=03 -O3" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -Wall -Mstandard -Mallocatable=03 -g" )
  
elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "PGI" )

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -Mstandard -Mallocatable=03 -fast" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -Mstandard -Mallocatable=03 -g" )
  
endif ()

add_subdirectory(src)

file(COPY Data DESTINATION ./ )
file(COPY input.dat DESTINATION ./ )

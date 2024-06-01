cmake_minimum_required(VERSION 3.16)

cmake_policy(SET CMP0069 NEW) # INTERPROCEDURAL_OPTIMIZATION
cmake_policy(SET CMP0083 NEW) # POSITION_INDEPENDENT_CODE

project(TDCS-e-2e Fortran)

set( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH} )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}" )
set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules )

if(NOT CMAKE_BUILD_TYPE)
  set( CMAKE_BUILD_TYPE "Release" )
  message(STATUS "No CMAKE_BUILD_TYPE selected, default to ${CMAKE_BUILD_TYPE}")
endif()

Option( TARGET "Generate target optimized binaries" ON )

if ( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" )

  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wimplicit-procedure -Wunused-parameter -std=f2018" )
  set( CMAKE_Fortran_FLAGS_FAST "-Ofast" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fcheck=all -ffpe-trap=zero,overflow,invalid -fbacktrace" )
  if( TARGET )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native" )
  endif()

elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel" )

  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -warn all -stand f18 -standard-semantics" )
  set( CMAKE_Fortran_FLAGS_FAST "-Ofast" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check all,nouninit -fpe-all=0 -traceback" )
  if( TARGET )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xHost" )
  endif()

else ()

  set( CMAKE_Fortran_FLAGS_RELEASE "-O3" )
  set( CMAKE_Fortran_FLAGS_DEBUG "-O0 -g" )

endif ()

find_package( OpenMP COMPONENTS Fortran)
option( ENABLE_OPENMP "Enable OpenMP" ${OpenMP_Fortran_FOUND} )

option( Profiling "Enable profiling with gprof" OFF )
if( Profiling )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -p" )
endif()

include( CheckIPOSupported )
check_ipo_supported( RESULT IPO_SUPPORTED OUTPUT IPO_MESSAGE LANGUAGES Fortran )
# LTO is broken with MINGW/binutils https://sourceware.org/bugzilla/show_bug.cgi?id=12762
if( MINGW AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
  set( IPO_SUPPORTED OFF )
endif()
if( IPO_SUPPORTED )
  set( CMAKE_INTERPROCEDURAL_OPTIMIZATION ON )
endif()

include( CheckPIESupported )
check_pie_supported( OUTPUT_VARIABLE PIE_MESSAGE LANGUAGES Fortran )
if( CMAKE_Fortran_LINK_PIE_SUPPORTED )
  set( CMAKE_POSITION_INDEPENDENT_CODE TRUE )
endif()

add_subdirectory(src)

file(COPY data DESTINATION ./ )
file(COPY input_example.dat DESTINATION ./ )

option( BUILD_TESTS "Build Tests" OFF)
if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()

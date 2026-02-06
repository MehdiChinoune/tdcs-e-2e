cmake_minimum_required(VERSION 3.22)

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

Option( TDCS_TARGET_NATIVE "Generate native target optimized binaries" ON )

if ( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" )

  if( CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 14.1 )
    set(_fstd 2023)
  else()
    set(_fstd 2018)
  endif()

  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wimplicit-procedure -Wunused-parameter -std=f${_fstd}" )
  set( CMAKE_Fortran_FLAGS_FAST "-Ofast" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fcheck=all -ffpe-trap=zero,overflow,invalid -fbacktrace" )
  if( TDCS_TARGET_NATIVE )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native" )
  endif()

elseif ( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "LLVMFlang" )

  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fimplicit-none-ext -std=f2018" )
  set( CMAKE_Fortran_FLAGS_FAST "-O3 -ffast-math -fstack-arrays" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -pedantic" )
  if( TDCS_TARGET_NATIVE )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native" )
  endif()

elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel" )

  if(WIN32)
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /warn:all /stand:f2023 -standard-semantics" )
    set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /check:all -fpe-all=0 -traceback" )
    if( TDCS_TARGET_NATIVE )
      set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /QxHost" )
    endif()
  else()
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -warn all -stand f2023 -standard-semantics" )
    set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check all -fpe-all=0 -traceback" )
    if( TDCS_TARGET_NATIVE )
      set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xHost" )
    endif()
  endif()
  set( CMAKE_Fortran_FLAGS_FAST "-Ofast" )

elseif ( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NVHPC" )

  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mallocatable=03 -Minform=warn -Mstandard" )
  set( CMAKE_Fortran_FLAGS_FAST "-fast" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Minfo=all -Mneginfo=all -Kieee" )
  #-Ktrap=fp" ) #for some reason it makes it read same line twice.
  if( TDCS_TARGET_NATIVE )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -tp host" )
  endif()

else ()

  set( CMAKE_Fortran_FLAGS_RELEASE "-O3" )
  set( CMAKE_Fortran_FLAGS_DEBUG "-O0 -g" )

endif ()

find_package( OpenMP COMPONENTS Fortran)
option( TDCS_OPENMP "Enable OpenMP" ${OpenMP_Fortran_FOUND} )

option( TDCS_PROFILE "Enable profiling with gprof" OFF )
if( TDCS_PROFILE )
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

add_custom_target( sync_data ALL
  COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different ${CMAKE_SOURCE_DIR}/data $<TARGET_FILE_DIR:main>/data
  COMMENT "Sync Data"
  )

file(COPY input_example.dat DESTINATION ./ )

option( TDCS_BUILD_TESTS "Build Tests" OFF)
option( TDCS_BUILD_SLOW_TESTS "Build Tests" OFF)
if(TDCS_BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()

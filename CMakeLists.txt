cmake_minimum_required(VERSION 3.7)

project(Main_new Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set( CMAKE_BUILD_TYPE "Release" )
  message(STATUS "No CMAKE_BUILD_TYPE selected, default to ${CMAKE_BUILD_TYPE}")
endif()

if ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" )

  if( "${CMAKE_Fortran_COMPILER_VERSION}" VERSION_LESS 6.1.0 )
    message("Unsupported Compiler you should use gfortran >= 6.1.0")
    return()
  endif()

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -Wall -Ofast" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -Wall -g -fcheck=all" )
  
elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "INTEL" )

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -warn all -Ofast -ipo -xHost" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -warn all -g" )
  
endif ()

add_subdirectory(src)

#set(srcs main.F90 constants.F90 trigo.F90 trigo_h.F90 special_functions_h.F90 special_functions.F90 utils_h.F90 utils.F90 fdcs_e2e_h.F90 fdcs_e2e.F90 )

#add_executable(main.exe ${srcs})

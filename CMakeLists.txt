cmake_minimum_required(VERSION 3.7)

project(Main_new Fortran)

# set ( CMAKE_Fortran_COMPILER "${CMAKE_Fortran_COMPILER}" )

if     ( "${CMAKE_Fortran_COMPILER}" MATCHES "gfortran" )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Ofast -flto -march=native" )
elseif ( "${CMAKE_Fortran_COMPILER}" MATCHES "ifort" )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -warn all -Ofast -ipo -march=native" )
endif ()

set(srcs main.F90 constants.F90 trigo.F90 trigo_h.F90 special_functions_h.F90 special_functions.F90 utils_h.F90 utils.F90 fdcs_e2e_h.F90 fdcs_e2e.F90 )

add_executable(main ${srcs})

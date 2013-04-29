project(lorentz)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH "/Users/phil/repos/cmake")
set(CMAKE_INSTALL_PREFIX /Users/phil/project/lorentz)
set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_Fortran_MODULE_DIRECTORY  ${CMAKE_INSTALL_PREFIX}/include_gfortran)
enable_language(Fortran)
set(CMAKE_Fortran_FLAGS "-fno-second-underscore -fno-range-check")

FILE(GLOB files_f90 "*.f90")

add_library(lorentz SHARED ${files_f90})

install(TARGETS lorentz LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX})



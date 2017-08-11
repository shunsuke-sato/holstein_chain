cmake_minimum_required(VERSION 2.8)

project(SPIN_BOSON Fortran)

#set(CMAKE_Fortran_COMPILER      "mpif90")
#set(CMAKE_Fortran_FLAGS   "-O2")

find_package(MPI REQUIRED)
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARIES})

#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})

set(TARGET_NAME   "holstain-chain")

set(BASIC_LIB "basic_lib")
set(MOD_LIB   "module_lib")
set(PBME_LIB  "pbme_lib")
set(FBTS_LIB  "fbts_lib")
set(CTEF_LIB  "ctef_lib")
set(PTEF_LIB  "ptef_lib")
set(GQME_LIB  "gqme_lib")
set(EXTERNAL_LIB "external_lib")

set(SRC_DIR "./src")
add_subdirectory(${SRC_DIR})


cmake_minimum_required(VERSION 2.8)

project(SPIN_BOSON Fortran)

set(CMAKE_Fortran_COMPILER      "mpif90")
set(CMAKE_Fortran_FLAGS   "-O2")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

set(TARGET_NAME   "holstain-chain")

set(SOURCES src/main.f90)
file(GLOB BASIC_SRC src/basic/*.f90)
file(GLOB PBME_SRC src/PBME/*.f90)
file(GLOB FBTS_SRC src/FBTS/*.f90)
file(GLOB PTEF_SRC src/PTEF/*.f90)
file(GLOB CTEF_SRC src/CTEF/*.f90)
file(GLOB GQME_SRC src/GQME/*.f90)
file(GLOB EXTERNAL_SRC src/external/*.f90)

list(APPEND SOURCES 
            ${BASIC_SRC}
            ${PBME_SRC}
            ${FBTS_SRC}
            ${PTEF_SRC}
            ${CTEF_SRC}
            ${GQME_SRC}
            ${EXTERNAL_SRC}
)


add_executable(${TARGET_NAME} ${SOURCES})

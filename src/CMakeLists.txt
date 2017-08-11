add_subdirectory("./basic")
add_subdirectory("./module")
add_subdirectory("./PBME")
add_subdirectory("./FBTS")
add_subdirectory("./CTEF")
add_subdirectory("./PTEF")
add_subdirectory("./GQME")
add_subdirectory("./external")

# Executable binary
add_executable(${TARGET_NAME} main.f90)
target_link_libraries(${TARGET_NAME} ${LIBRARIES} ${MPI_Fortran_LIBRARIES})
add_dependencies(${TARGET_NAME} ${LIBRARIES})

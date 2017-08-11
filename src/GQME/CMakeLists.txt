set(TARGET_NAME ${GQME_LIB})
set(SOURCES
GQME_dynamics.f90
GQME_kernel.f90
GQME_kernel_read.f90
refine_GQME_kernel_full.f90
refine_GQME_kernel_K1K3.f90
)

add_library(${TARGET_NAME} STATIC ${SOURCES})
add_dependencies(${TARGET_NAME} ${BASIC_LIB}
                                ${MOD_LIB})
list(APPEND LIBRARIES ${TARGET_NAME})
set(LIBRARIES ${LIBRARIES} PARENT_SCOPE)



set(SOURCES
FBTS_allocation.f90
FBTS_dt_evolve.f90
FBTS_dynamics.f90
FBTS_initial_condition.f90
FBTS_kernel.f90
)

add_library(${FBTS_LIB} STATIC ${SOURCES})
add_dependencies(${FBTS_LIB} ${BASIC_LIB}
                                ${MOD_LIB})
list(APPEND LIBRARIES ${FBTS_LIB})
set(LIBRARIES ${LIBRARIES} PARENT_SCOPE)



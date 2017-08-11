set(SOURCES
PBME_allocation.f90
PBME_dt_evolve.f90
PBME_dynamics.f90
PBME_initial_condition.f90
PBME_kernel.f90
)

add_library(${PBME_LIB} STATIC ${SOURCES})
add_dependencies(${PBME_LIB} ${BASIC_LIB})
list(APPEND LIBRARIES ${PBME_LIB})
set(LIBRARIES ${LIBRARIES} PARENT_SCOPE)


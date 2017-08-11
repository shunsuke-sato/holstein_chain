set(SOURCES
PTEF_allocation.f90
PTEF_calc_matrix.f90
PTEF_dynamics.f90
PTEF_initial_condition.f90
)

add_library(${PTEF_LIB} STATIC ${SOURCES})
add_dependencies(${PTEF_LIB} ${BASIC_LIB})
list(APPEND LIBRARIES ${PTEF_LIB})
set(LIBRARIES ${LIBRARIES} PARENT_SCOPE)



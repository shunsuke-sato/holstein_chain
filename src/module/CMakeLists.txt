set(SOURCES
communication.f90
parallel.f90
math_module.f90
)

add_library(${MOD_LIB} STATIC ${SOURCES})
list(APPEND LIBRARIES ${MOD_LIB})
set(LIBRARIES ${LIBRARIES} PARENT_SCOPE)


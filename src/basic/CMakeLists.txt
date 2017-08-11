set(SOURCES
      allocate_GQME_kernel.f90
      allocation_of_global_arrays.f90
      calc_energy.f90
      calc_force_HO.f90
      correlated_gaussian_random_number.f90
      dt_evolve_elec.f90
      dt_evolve_elec_pair.f90
      err_finalize.f90
      evaluate_GQME_kernel_full.f90
      evaluate_GQME_kernel_K1K3.f90
      gaussian_random_number.f90
      global_variables.f90
      hpsi.f90
      initialize_mpi.f90
      initialize_random_number_generator.f90
      multi_traject_Ehrenfest.f90
      ranlux_double.f90
      set_initial_conditions.f90
      set_initial_conditions_elec.f90
      set_initial_conditions_ph.f90
      set_model_parameters.f90
      set_thermal_ph_dist.f90
)

add_library(${BASIC_LIB} STATIC ${SOURCES})
add_dependencies(${BASIC_LIB} ${MOD_LIB})
list(APPEND LIBRARIES ${BASIC_LIB})
set(LIBRARIES ${LIBRARIES} PARENT_SCOPE)

set(RegistrySources
        register_memory_pool.cpp
        register_refinement_algorithms.cpp
        register_preprocessing_algorithms.cpp
        register_initial_partitioning_algorithms.cpp
        register_policies.cpp
        register_coarsening_algorithms.cpp
        )


foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
    target_sources(${modtarget} PRIVATE ${RegistrySources})
endforeach()
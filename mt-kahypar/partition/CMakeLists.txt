add_subdirectory(preprocessing)
add_subdirectory(coarsening)
add_subdirectory(refinement)
add_subdirectory(initial_partitioning)
add_subdirectory(mapping)
add_subdirectory(registries)

set(PartitionSources
        partitioner.cpp
        partitioner_facade.cpp
        multilevel.cpp
        context.cpp
        context_enum_classes.cpp
        conversion.cpp
        metrics.cpp
        recursive_bipartitioning.cpp
        )

set(DeepMultilevelSources
        deep_multilevel.cpp
        )

target_sources(MtKaHyPar-Sources INTERFACE ${PartitionSources})

if (NOT KAHYPAR_MINIMAL_COMPILATION)
  target_sources(MtKaHyPar-Sources INTERFACE ${DeepMultilevelSources})
endif()

set(PartitionToolsSources
        context.cpp
        context_enum_classes.cpp
        conversion.cpp
        metrics.cpp
        )

target_sources(MtKaHyPar-ToolsSources INTERFACE ${PartitionToolsSources})

add_subdirectory(preprocessing/)
add_subdirectory(coarsening)
add_subdirectory(refinement)
add_subdirectory(initial_partitioning)
add_subdirectory(registries)

set(PartitionSources
        partitioner.cpp
        multilevel.cpp
        context.cpp
        context_enum_classes.cpp
        metrics.cpp
        recursive_bipartitioning.cpp
        deep_multilevel.cpp
        )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
    target_sources(${modtarget} PRIVATE ${PartitionSources})
endforeach()

set(ContextWrapperSources
        context.cpp
        context_enum_classes.cpp
        )
target_sources(MtKaHyParWrapper PRIVATE ${ContextWrapperSources})

set(ContextOneMtKaHyParSources
        context.cpp
        metrics.cpp
        context_enum_classes.cpp
        partitioner_facade.cpp
        )
target_sources(OneMtKaHyPar PRIVATE ${ContextOneMtKaHyParSources})

add_subdirectory(preprocessing/)
add_subdirectory(coarsening)
add_subdirectory(refinement)

set(PartitionSources
        partitioner.cpp
        multilevel.cpp
        context.cpp
        context_enum_classes.cpp
        )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
    target_sources(${modtarget} PRIVATE ${PartitionSources})
endforeach()
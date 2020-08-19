add_subdirectory(preprocessing/)
add_subdirectory(coarsening)
add_subdirectory(refinement)

set(PartitionSources

        )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
    target_sources(${modtarget} PRIVATE ${PartitionSources})
endforeach()
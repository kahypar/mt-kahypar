set(FixedVertexSources
      fixed_vertex_removal.cpp)

foreach(modtarget IN LISTS PARTITIONING_SUITE_TARGETS)
    target_sources(${modtarget} PRIVATE ${FixedVertexSources})
endforeach()
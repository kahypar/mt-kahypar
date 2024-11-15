#################################################################
## Minimal tools with single source file                       ##
#################################################################

# meta target for single executable tools
add_library(MtKaHyPar-ToolsMin INTERFACE)
target_link_libraries(MtKaHyPar-ToolsMin INTERFACE MtKaHyPar-Include MtKaHyPar-BuildFlags MtKaHyPar-WarningFlags)

add_executable(SnapToHgr snap_graph_to_hgr.cc)
target_link_libraries(SnapToHgr MtKaHyPar-ToolsMin)

add_executable(SnapToMetis snap_to_metis_converter.cc)
target_link_libraries(SnapToMetis MtKaHyPar-ToolsMin)

add_executable(BenchShuffle bench_deterministic_shuffling.cpp)
target_link_libraries(BenchShuffle MtKaHyPar-ToolsMin)

add_executable(MtxToGraph mtx_to_graph.cc)
target_link_libraries(MtxToGraph MtKaHyPar-BuildFlags MtKaHyPar-WarningFlags)


#################################################################
## Tools with more dependencies                                ##
#################################################################

# meta target for tools with more dependencies
add_library(MtKaHyPar-BuildTools OBJECT "")
target_link_libraries(MtKaHyPar-BuildTools PRIVATE MtKaHyPar-ToolsSources ${HWLOCK_LIBRARY})
target_link_libraries(MtKaHyPar-BuildTools PUBLIC MtKaHyPar-Include MtKaHyPar-BuildFlags MtKaHyPar-WarningFlags)

add_executable(GraphToHgr graph_to_hgr.cc)
target_link_libraries(GraphToHgr MtKaHyPar-BuildTools)

add_executable(HgrToGraph hgr_to_graph.cc)
target_link_libraries(HgrToGraph MtKaHyPar-BuildTools)

add_executable(HgrToParkway hgr_to_parkway_converter.cc)
target_link_libraries(HgrToParkway MtKaHyPar-BuildTools)

add_executable(HgrToZoltan hgr_to_zoltan_converter.cc)
target_link_libraries(HgrToZoltan MtKaHyPar-BuildTools)

add_executable(MetisToScotch metis_to_scotch_converter.cc)
target_link_libraries(MetisToScotch MtKaHyPar-BuildTools)

add_executable(HypergraphStats hypergraph_stats.cc)
target_link_libraries(HypergraphStats MtKaHyPar-BuildTools)

add_executable(OneToOneMapping one_to_one_mapping.cc)
target_link_libraries(OneToOneMapping MtKaHyPar-BuildTools)

add_executable(VerifyTargetGraphPartition verify_process_mapping_partition.cc)
target_link_libraries(VerifyTargetGraphPartition MtKaHyPar-BuildTools)

add_executable(GridGraphGenerator grid_graph_generator.cc)
target_link_libraries(GridGraphGenerator MtKaHyPar-BuildTools)

add_executable(HierarchicalTargetGraphGenerator hierarchical_process_graph_generator.cc)
target_link_libraries(HierarchicalTargetGraphGenerator MtKaHyPar-BuildTools)

add_executable(EvaluateBipart evaluate_bipart_partition.cc)
target_link_libraries(EvaluateBipart MtKaHyPar-BuildTools)

add_executable(EvaluatePartition evaluate_hmetis_partition.cc)
target_link_libraries(EvaluatePartition MtKaHyPar-BuildTools)

add_executable(VerifyPartition verify_partition.cc)
target_link_libraries(VerifyPartition MtKaHyPar-BuildTools)

add_executable(FixedVertexFileGenerator fixed_vertex_file_generator.cc)
target_link_libraries(FixedVertexFileGenerator MtKaHyPar-BuildTools)

add_executable(SnapGraphToHgr snap_graph_to_hgr.cc)
target_link_libraries(SnapGraphToHgr ${Boost_LIBRARIES})
set_property(TARGET SnapGraphToHgr PROPERTY CXX_STANDARD 17)
set_property(TARGET SnapGraphToHgr PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(GraphToHgr graph_to_hgr.cc)
target_link_libraries(GraphToHgr ${Boost_LIBRARIES})
set_property(TARGET GraphToHgr PROPERTY CXX_STANDARD 17)
set_property(TARGET GraphToHgr PROPERTY CXX_STANDARD_REQUIRED ON)


add_executable(HgrToParkway hgr_to_parkway_converter.cc)
target_link_libraries(HgrToParkway ${Boost_LIBRARIES})
set_property(TARGET HgrToParkway PROPERTY CXX_STANDARD 17)
set_property(TARGET HgrToParkway PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(HgrToZoltan hgr_to_zoltan_converter.cc)
target_link_libraries(HgrToZoltan ${Boost_LIBRARIES})
set_property(TARGET HgrToZoltan PROPERTY CXX_STANDARD 17)
set_property(TARGET HgrToZoltan PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(HypergraphStats hypergraph_stats.cc)
target_link_libraries(HypergraphStats ${Boost_LIBRARIES})
set_property(TARGET HypergraphStats PROPERTY CXX_STANDARD 17)
set_property(TARGET HypergraphStats PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(MetisToScotch metis_to_scotch_converter.cc)
target_link_libraries(MetisToScotch ${Boost_LIBRARIES})
set_property(TARGET MetisToScotch PROPERTY CXX_STANDARD 17)
set_property(TARGET MetisToScotch PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(SnapToMetis snap_to_metis_converter.cc)
target_link_libraries(SnapToMetis ${Boost_LIBRARIES})
set_property(TARGET SnapToMetis PROPERTY CXX_STANDARD 17)
set_property(TARGET SnapToMetis PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(EvaluateBipart evaluate_bipart_partition.cc)
target_link_libraries(EvaluateBipart ${Boost_LIBRARIES})
set_property(TARGET EvaluateBipart PROPERTY CXX_STANDARD 17)
set_property(TARGET EvaluateBipart PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(EvaluatePartition evaluate_hmetis_partition.cc)
target_link_libraries(EvaluatePartition ${Boost_LIBRARIES})
set_property(TARGET EvaluatePartition PROPERTY CXX_STANDARD 17)
set_property(TARGET EvaluatePartition PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(VerifyPartition verify_partition.cc)
target_link_libraries(VerifyPartition ${Boost_LIBRARIES})
set_property(TARGET VerifyPartition PROPERTY CXX_STANDARD 17)
set_property(TARGET VerifyPartition PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(BenchShuffle bench_deterministic_shuffling.cpp bench_deterministic_shuffling.cpp)
set_property(TARGET BenchShuffle PROPERTY CXX_STANDARD 17)
set_property(TARGET BenchShuffle PROPERTY CXX_STANDARD_REQUIRED ON)

set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} EvaluateBipart EvaluatePartition VerifyPartition HgrToZoltan HypergraphStats MetisToScotch SnapToMetis GraphToHgr HgrToParkway SnapGraphToHgr PARENT_SCOPE)

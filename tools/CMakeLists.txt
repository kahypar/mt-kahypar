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

add_executable(EvaluateBipart evaluate_bipart_partition.cc)
target_link_libraries(EvaluateBipart ${Boost_LIBRARIES})
set_property(TARGET EvaluateBipart PROPERTY CXX_STANDARD 17)
set_property(TARGET EvaluateBipart PROPERTY CXX_STANDARD_REQUIRED ON)

set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} EvaluateBipart HgrToZoltan GraphToHgr HgrToParkway SnapGraphToHgr PARENT_SCOPE)
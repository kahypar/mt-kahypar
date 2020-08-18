add_executable(SnapGraphToHgr snap_graph_to_hgr.cc snap_graph_to_hgr.cc)
set_property(TARGET SnapGraphToHgr PROPERTY CXX_STANDARD 14)
set_property(TARGET SnapGraphToHgr PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(SnapGraphToHgr ${Boost_LIBRARIES})

add_executable(GraphToHgr graph_to_hgr.cc graph_to_hgr.cc)
set_property(TARGET GraphToHgr PROPERTY CXX_STANDARD 14)
set_property(TARGET GraphToHgr PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(GraphToHgr ${Boost_LIBRARIES})


add_executable(HgrToParkway hgr_to_parkway_converter.cc hgr_to_parkway_converter.cc)
set_property(TARGET HgrToParkway PROPERTY CXX_STANDARD 14)
set_property(TARGET HgrToParkway PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(HgrToParkway ${Boost_LIBRARIES})

add_executable(HgrToZoltan hgr_to_zoltan_converter.cc hgr_to_zoltan_converter.cc)
set_property(TARGET HgrToZoltan PROPERTY CXX_STANDARD 14)
set_property(TARGET HgrToZoltan PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(HgrToZoltan ${Boost_LIBRARIES})

add_executable(SnapGraphToHgr snap_graph_to_hgr.cc snap_graph_to_hgr.cc)
set_property(TARGET SnapGraphToHgr PROPERTY CXX_STANDARD 14)
set_property(TARGET SnapGraphToHgr PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(SnapGraphToHgr ${Boost_LIBRARIES})
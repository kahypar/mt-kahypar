################################################################################
# python/CMakeLists.txt
#
# Copyright (c) 2019 Timo Bingmann
#
# All rights reserved. Published under the MIT License in the LICENSE file.
################################################################################

add_subdirectory(pybind11)
include_directories(${PROJECT_SOURCE_DIR})
pybind11_add_module(mtkahypar_python module.cpp)

target_link_libraries(mtkahypar_python PRIVATE TBB::tbb TBB::tbbmalloc_proxy)

if(NOT MT_KAHYPAR_DISABLE_BOOST)
    target_link_libraries(mtkahypar_python PRIVATE Boost::program_options)
endif()

SET_TARGET_PROPERTIES(mtkahypar_python PROPERTIES COMPILE_FLAGS "-fvisibility=hidden")
target_compile_definitions(mtkahypar_python PUBLIC MT_KAHYPAR_LIBRARY_MODE)
target_compile_definitions(mtkahypar_python PUBLIC KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES)
target_compile_definitions(mtkahypar_python PUBLIC KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES)
target_compile_definitions(mtkahypar_python PUBLIC KAHYPAR_ENABLE_SOED_METRIC)
target_compile_definitions(mtkahypar_python PUBLIC KAHYPAR_ENABLE_STEINER_TREE_METRIC)

set(PARTITIONING_SUITE_TARGETS ${PARTITIONING_SUITE_TARGETS} mtkahypar_python PARENT_SCOPE)

# rename mtkahypar_python target output to mtkahypar
set_target_properties(mtkahypar_python PROPERTIES OUTPUT_NAME mtkahypar)

################################################################################

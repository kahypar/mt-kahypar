################################################################################
# python/CMakeLists.txt
#
# Copyright (c) 2019 Timo Bingmann
#
# All rights reserved. Published under the MIT License in the LICENSE file.
################################################################################

add_subdirectory(pybind11)
include_directories(${PROJECT_SOURCE_DIR})
pybind11_add_module(mtkahyparhgp_python module_hgp.cpp)
target_link_libraries(mtkahyparhgp_python PRIVATE ${Boost_LIBRARIES})
SET_TARGET_PROPERTIES(mtkahyparhgp_python PROPERTIES COMPILE_FLAGS "-fvisibility=hidden")
target_compile_definitions(mtkahyparhgp_python PUBLIC MT_KAHYPAR_LIBRARY_MODE)

pybind11_add_module(mtkahypargp_python module_gp.cpp)
target_link_libraries(mtkahypargp_python PRIVATE ${Boost_LIBRARIES})
SET_TARGET_PROPERTIES(mtkahypargp_python PROPERTIES COMPILE_FLAGS "-fvisibility=hidden")
target_compile_definitions(mtkahypargp_python PUBLIC USE_GRAPH_PARTITIONER)
target_compile_definitions(mtkahypargp_python PUBLIC MT_KAHYPAR_LIBRARY_MODE)

set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} mtkahyparhgp_python mtkahypargp_python PARENT_SCOPE)

# rename mtkahypar_python target output to mtkahypar
set_target_properties(mtkahyparhgp_python PROPERTIES OUTPUT_NAME mtkahyparhgp)
set_target_properties(mtkahypargp_python PROPERTIES OUTPUT_NAME mtkahypargp)

################################################################################

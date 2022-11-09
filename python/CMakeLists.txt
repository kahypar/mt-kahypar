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
target_link_libraries(mtkahypar_python PRIVATE ${Boost_LIBRARIES})

# rename mtkahypar_python target output to mtkahypar
set_target_properties(mtkahypar_python PROPERTIES OUTPUT_NAME mtkahypar)

################################################################################

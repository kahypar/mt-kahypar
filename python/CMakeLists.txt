add_subdirectory(pybind11)
# include_directories(${PROJECT_SOURCE_DIR})
pybind11_add_module(mtkahypar_python module.cpp)

target_link_libraries(mtkahypar_python PRIVATE MtKaHyPar-LibraryBuildSources)

# TODO
SET_TARGET_PROPERTIES(mtkahypar_python PROPERTIES COMPILE_FLAGS "-fvisibility=hidden")

# rename mtkahypar_python target output to mtkahypar
set_target_properties(mtkahypar_python PROPERTIES OUTPUT_NAME mtkahypar)


FetchContent_Declare(
    pybind11 SYSTEM EXCLUDE_FROM_ALL
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        ${KAHYPAR_PYBIND11_VERSION}
)
FetchContent_MakeAvailable(pybind11)
pybind11_add_module(mtkahypar_python module.cpp)

target_link_libraries(mtkahypar_python PRIVATE MtKaHyPar-LibraryBuildSources)
set_target_properties(mtkahypar_python PROPERTIES COMPILE_FLAGS "-fvisibility=hidden")

# rename mtkahypar_python target output to mtkahypar
set_target_properties(mtkahypar_python PROPERTIES OUTPUT_NAME mtkahypar)

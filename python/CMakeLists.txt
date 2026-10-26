FetchContent_Declare(
    pybind11 SYSTEM EXCLUDE_FROM_ALL
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        ${KAHYPAR_PYBIND11_VERSION}
)
set(PYBIND11_FINDPYTHON ON)
FetchContent_MakeAvailable(pybind11)

pybind11_add_module(mtkahypar_python module.cpp)
hide_symbols(mtkahypar_python)
copy_dependencies(mtkahypar_python)

target_link_libraries(mtkahypar_python PRIVATE MtKaHyPar-LibraryBuildSources)

# rename mtkahypar_python target output to mtkahypar
set_target_properties(mtkahypar_python PROPERTIES
    INTERPROCEDURAL_OPTIMIZATION ON
    OUTPUT_NAME mtkahypar)

add_subdirectory(tests)

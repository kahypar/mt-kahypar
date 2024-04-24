# file(GLOB SOLVER solvers/julia_gevp.cpp)
file(GLOB SOLVER solvers/slepc_gevp.cpp)
set(AlgebraicWrappers
    algebraic_wrappers/vector.cpp
    algebraic_wrappers/operator.cpp
)

foreach(modtarget IN LISTS PARTITIONING_SUITE_TARGETS)
    target_sources(${modtarget} PRIVATE ${Solver})
    target_sources(${modtarget} PRIVATE ${AlgebraicWrappers})
    target_precompile_headers(${modtarget} PRIVATE datatypes.h)
endforeach()

# set_source_files_properties(${Solver} PROPERTIES COMPILE_OPTIONS -std=c++20)

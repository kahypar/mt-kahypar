add_subdirectory(partition)
add_subdirectory(parallel)
add_subdirectory(datastructures)
add_subdirectory(io)
add_subdirectory(utils)

target_precompile_headers(MtKaHyPar-Sources INTERFACE definitions.h)
target_precompile_headers(MtKaHyPar-Sources INTERFACE partition/context.h)

add_executable(MtKaHyPar mt_kahypar.cc)
target_link_libraries(MtKaHyPar MtKaHyPar-BuildSources)


# additional installation target since installation-time renaming does not seem to be possible
add_executable(MtKaHyPar-CLI mt_kahypar.cc)
target_link_libraries(MtKaHyPar-CLI MtKaHyPar-BuildSources)

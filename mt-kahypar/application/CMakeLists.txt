<<<<<<< HEAD
add_mtkahypar_executable(MtKaHyPar mt_kahypar.cc)
=======
add_executable(MtKaHyPar mt_kahypar.cc)
add_executable(MtKaHyParC compression.cc)
>>>>>>> origin/compression
target_link_libraries(MtKaHyPar MtKaHyPar-BuildSources)
target_link_libraries(MtKaHyParC MtKaHyPar-BuildSources)


# additional installation target since installation-time renaming does not seem to be possible
add_mtkahypar_executable(MtKaHyPar-CLI mt_kahypar.cc)
target_link_libraries(MtKaHyPar-CLI MtKaHyPar-BuildSources)

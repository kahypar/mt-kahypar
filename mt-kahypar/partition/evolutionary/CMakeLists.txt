set(EvolutionarySources
        action.h
        combine.h
        diversifier.h
        edge_frequency.h
        evo_logs.h
        individual.h
        individual.cpp
        mutate.h
        population.h
        population.cpp
        strategy_picker.h
        strategy_picker.cpp
        stopping_criteria.h
        stopping_criteria.cpp)

foreach(modtarget IN LISTS PARTITIONING_SUITE_TARGETS)
    target_sources(${modtarget} PRIVATE ${EvolutionarySources})
endforeach()

# Ensure the application build target also gets these sources
target_sources(MtKaHyPar-BuildSources INTERFACE ${EvolutionarySources})
target_sources(MtKaHyPar-ToolsSources INTERFACE ${EvolutionarySources})

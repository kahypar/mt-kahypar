set(EvolutionarySources
        action.h
        combine.h
        diversifier.h
        edge_frequency.h
        individual.h
        individual.cpp
        mutate.h
        population.h
        population.cpp
        probability_tables.h)

foreach(modtarget IN LISTS PARTITIONING_SUITE_TARGETS)
    target_sources(${modtarget} PRIVATE ${EvolutionarySources})
endforeach()

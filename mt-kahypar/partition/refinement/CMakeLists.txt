set(RefinementSources
        fm/fm_commons.cpp
        fm/multitry_kway_fm.cpp
        fm/localized_kway_fm_core.cpp
        fm/global_rollback.cpp
        gains/km1/km1_gain_cache.cpp
        gains/cut/cut_gain_cache.cpp
        gains/soed/soed_gain_cache.cpp
        gains/cut_for_graphs/cut_gain_cache_for_graphs.cpp
        fm/sequential_twoway_fm_refiner.cpp
        label_propagation/label_propagation_refiner.cpp
        jet/jet_refiner.cpp
        rebalancing/rebalancer.cpp
        rebalancing/jet_rebalancer.cpp
        rebalancing/kaminpar_rebalancer.cpp
        rebalancing/call_kaminpar.cpp
        deterministic/deterministic_label_propagation.cpp
        flows/refiner_adapter.cpp
        flows/problem_construction.cpp
        flows/scheduler.cpp
        flows/quotient_graph.cpp
        flows/flow_refiner.cpp
        flows/sequential_construction.cpp
        flows/parallel_construction.cpp
        flows/flow_hypergraph_builder.cpp
        )

foreach(modtarget IN LISTS PARTITIONING_SUITE_TARGETS)
    target_sources(${modtarget} PRIVATE ${RefinementSources})
endforeach()
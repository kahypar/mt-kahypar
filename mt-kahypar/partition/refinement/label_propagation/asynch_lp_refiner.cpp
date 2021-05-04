//
// Created by mlaupichler on 04.05.21.
//

#include "asynch_lp_refiner.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

    template<template <typename> class GainPolicy>
    bool AsynchLPRefiner<GainPolicy>::refineImpl(PartitionedHypergraph &hypergraph,
                                                 const parallel::scalable_vector <mt_kahypar::HypernodeID> &refinement_nodes,
                                                 kahypar::Metrics &best_metrics, double) {

        ASSERT(!refinement_nodes.empty(), "AsynchLPRefiner will not work without given seed refinement nodes. Cannot be used "
                                          "solely for rebalancing or for global refinement!");
        _active_nodes = refinement_nodes;

        _gain.reset();
        _next_active.reset();

        // Perform Label Propagation
        labelPropagation(hypergraph);

        // Update global part weight and sizes
        best_metrics.imbalance = metrics::imbalance(hypergraph, _context);

        // Update metrics statistics
        HyperedgeWeight current_metric = best_metrics.getMetric(
                kahypar::Mode::direct_kway, _context.partition.objective);
        Gain delta = _gain.delta();
        ASSERT(delta <= 0, "LP refiner worsen solution quality");

        HEAVY_REFINEMENT_ASSERT(hypergraph.checkTrackedPartitionInformation());
        HEAVY_REFINEMENT_ASSERT(current_metric + delta ==
                                metrics::objective(hypergraph, _context.partition.objective,
                                                   !_context.refinement.label_propagation.execute_sequential),
                                V(current_metric) << V(delta) <<
                                                  V(metrics::objective(hypergraph, _context.partition.objective,
                                                                       _context.refinement.label_propagation.execute_sequential)));

        best_metrics.updateMetric(current_metric + delta, kahypar::Mode::direct_kway, _context.partition.objective);
        utils::Stats::instance().update_stat("lp_improvement", std::abs(delta));
        return delta < 0;

    }

    template <template <typename> class GainPolicy>
    void AsynchLPRefiner<GainPolicy>::labelPropagation(PartitionedHypergraph& hypergraph) {
        NextActiveNodes next_active_nodes;
        for (size_t i = 0; i < _context.refinement.label_propagation.maximum_iterations; ++i) {
            DBG << "Starting Label Propagation Round" << i;

            utils::Timer::instance().start_timer(
                    "lp_round_" + std::to_string(i), "Label Propagation Round " + std::to_string(i), true);

            if ( !_active_nodes.empty() ) {
                labelPropagationRound(hypergraph, next_active_nodes);
            }

            if ( _context.refinement.label_propagation.execute_sequential ) {
                _active_nodes = next_active_nodes.copy_sequential();
                next_active_nodes.clear_sequential();
            } else {
                _active_nodes = next_active_nodes.copy_parallel();
                next_active_nodes.clear_parallel();
            }
            utils::Timer::instance().stop_timer("lp_round_" + std::to_string(i));

            if ( _active_nodes.empty() ) {
                break;
            }
        }
    }

    template <template <typename> class GainPolicy>
    bool AsynchLPRefiner<GainPolicy>::labelPropagationRound(
            PartitionedHypergraph& hypergraph,
            NextActiveNodes& next_active_nodes) {

        _visited_he.reset();
        _next_active.reset();
        // This function is passed as lambda to the changeNodePart function and used
        // to calculate the "real" delta of a move (in terms of the used objective function).
        auto objective_delta = [&](const HyperedgeID he,
                                   const HyperedgeWeight edge_weight,
                                   const HypernodeID edge_size,
                                   const HypernodeID pin_count_in_from_part_after,
                                   const HypernodeID pin_count_in_to_part_after) {
            _gain.computeDeltaForHyperedge(he, edge_weight, edge_size,
                                           pin_count_in_from_part_after, pin_count_in_to_part_after);
        };

        // Shuffle Vector
        bool converged = true;
        if ( _context.refinement.label_propagation.execute_sequential ) {
            utils::Randomize::instance().shuffleVector(
                    _active_nodes, 0UL, _active_nodes.size(), sched_getcpu());

            for ( size_t j = 0; j < _active_nodes.size(); ++j ) {
                const HypernodeID hn = _active_nodes[j];
                if ( ! moveVertex(hypergraph, hn, next_active_nodes, objective_delta) ) {
                    converged = false;
                }
            }
        } else {
            utils::Randomize::instance().parallelShuffleVector(
                    _active_nodes, 0UL, _active_nodes.size());

            tbb::parallel_for(0UL, _active_nodes.size(), [&](const size_t& j) {
                const HypernodeID hn = _active_nodes[j];
                if ( ! moveVertex(hypergraph, hn, next_active_nodes, objective_delta) ) {
                    converged = false;
                }
            });
        }

        HEAVY_REFINEMENT_ASSERT(hypergraph.checkTrackedPartitionInformation());
        return converged;
    }

    // explicitly instantiate so the compiler can generate them when compiling this cpp file
    template class AsynchLPRefiner<Km1Policy>;
    template class AsynchLPRefiner<CutPolicy>;

} // namespace mt_kahypar



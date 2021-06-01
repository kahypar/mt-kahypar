//
// Created by mlaupichler on 04.05.21.
//

#include <tbb/parallel_sort.h>
#include "async_lp_refiner.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

    template<template <typename> class LocalGainPolicy>
    bool AsyncLPRefiner<LocalGainPolicy>::refineImpl(PartitionedHypergraph &hypergraph,
                                                     const parallel::scalable_vector <mt_kahypar::HypernodeID> &refinement_nodes,
                                                     metrics::ThreadSafeMetrics &best_metrics,
                                                     double,
                                                     ds::StreamingVector<HypernodeID>& moved_nodes) {
        ASSERT(_contraction_group_id != ds::invalidGroupID, "ContractionGroupID (Owner-ID) for locking is invalid.");
        ASSERT(!refinement_nodes.empty(), "AsyncLPRefiner will not work without given seed refinement nodes. Cannot be used "
                                          "solely for rebalancing or for global refinement!");
        ASSERT(std::all_of(refinement_nodes.begin(),refinement_nodes.end(),[&](const HypernodeID& hn) {return hypergraph.nodeIsEnabled(hn);})
               && "Not all given seed nodes are enabled!");

        _active_nodes.assign(refinement_nodes.begin(), refinement_nodes.end());

        _gain.reset();

        // Perform Label Propagation
        labelPropagation(hypergraph, moved_nodes);

        // Update global part weight and sizes
        double imbalance;
        do {
            imbalance = metrics::imbalance(hypergraph, _context);
        } while (!best_metrics.update_imbalance_strong(imbalance));

        // Update metrics statistics
        HyperedgeWeight current_metric = best_metrics.getMetric(
                kahypar::Mode::direct_kway, _context.partition.objective);
        Gain delta = _gain.delta();
        ASSERT(delta <= 0, "LP refiner worsen solution quality");

//        HEAVY_REFINEMENT_ASSERT(hypergraph.checkTrackedPartitionInformation());
//        HEAVY_REFINEMENT_ASSERT(current_metric + delta == metrics::objective(hypergraph, _context.partition.objective, false),
//                                V(current_metric) << V(delta) <<
//                                                  V(metrics::objective(hypergraph, _context.partition.objective, false)));

        best_metrics.fetch_add(delta, kahypar::Mode::direct_kway, _context.partition.objective);
        utils::Stats::instance().update_stat("lp_improvement", std::abs(delta));
        return delta < 0;

    }

    template <template <typename> class LocalGainPolicy>
    void AsyncLPRefiner<LocalGainPolicy>::labelPropagation(PartitionedHypergraph &hypergraph, MovedNodes &moved_nodes) {
        NextActiveNodes next_active_nodes;
        VisitedEdges visited_edges;
        for (size_t i = 0; i < _context.refinement.label_propagation.maximum_iterations; ++i) {
            DBG << "Starting Label Propagation Round" << i;

            utils::Timer::instance().start_timer(
                    "lp_round_" + std::to_string(i), "Label Propagation Round " + std::to_string(i), true);

            if ( !_active_nodes.empty() ) {
                labelPropagationRound(hypergraph, next_active_nodes, visited_edges, moved_nodes);
            }

            _active_nodes = next_active_nodes;

            auto has_node_duplicates = [&]() -> bool {
                std::sort(next_active_nodes.begin(),next_active_nodes.end());
                return std::adjacent_find(next_active_nodes.begin(), next_active_nodes.end()) != next_active_nodes.end();
            };
            ASSERT(! has_node_duplicates());

            auto has_edge_duplicates = [&]() -> bool {
                std::sort(visited_edges.begin(), visited_edges.end());
                return std::adjacent_find(visited_edges.begin(), visited_edges.end()) != visited_edges.end();
            };
            ASSERT(! has_edge_duplicates());

            // Linear reset of anti-duplicator flags that were set in the last round
            for (auto hn : next_active_nodes) {
                bool reset = _next_active.compare_and_set_to_false(hn);
                ASSERT(reset);
            }
            for (auto he : visited_edges) {
                bool reset = _visited_he.compare_and_set_to_false(he);
                ASSERT(reset);
            }

            next_active_nodes.clear();
            visited_edges.clear();

            utils::Timer::instance().stop_timer("lp_round_" + std::to_string(i));

            if ( _active_nodes.empty() ) {
                break;
            }
        }
    }

    template <template <typename> class LocalGainPolicy>
    bool AsyncLPRefiner<LocalGainPolicy>::labelPropagationRound(PartitionedHypergraph &hypergraph,
                                                                NextActiveNodes &next_active_nodes,
                                                                VisitedEdges &visited_edges, MovedNodes &moved_nodes) {

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
        std::shuffle(_active_nodes.begin(),_active_nodes.end(), _rng);

        bool converged = true;
        parallel::scalable_vector<HypernodeID> retry_nodes;
        for ( size_t j = 0; j < _active_nodes.size(); ++j ) {
            const HypernodeID hn = _active_nodes[j];
            bool acquired = _lock_manager->tryToAcquireLock(hn, _contraction_group_id);
            if (acquired) {
                if ( moveVertex(hypergraph, hn, next_active_nodes, visited_edges, objective_delta) ) {
                    moved_nodes.stream(hn);
                } else {
                    converged = false;
                }
                _lock_manager->strongReleaseLock(hn, _contraction_group_id);
            } else {
                retry_nodes.push_back(hn);
            }
        }

        // Retry acquiring lock and moving exactly once
        // todo mlaupichler: This is arbitrary, perhaps add CL option to set how often we retry here
        for ( size_t j = 0; j < retry_nodes.size(); ++j ) {
            const HypernodeID hn = retry_nodes[j];
            bool acquired = _lock_manager->tryToAcquireLock(hn, _contraction_group_id);
            if (acquired) {
                if ( moveVertex(hypergraph, hn, next_active_nodes, visited_edges, objective_delta) ) {
                    moved_nodes.stream(hn);
                } else {
                    converged = false;
                }
                _lock_manager->strongReleaseLock(hn, _contraction_group_id);
            }
        }

//        HEAVY_REFINEMENT_ASSERT(hypergraph.checkTrackedPartitionInformation());
        return converged;
    }

    template <template <typename> class LocalGainPolicy>
    void AsyncLPRefiner<LocalGainPolicy>::resetForGroup(ds::ContractionGroupID groupID) {
        _contraction_group_id = groupID;
    }

    // explicitly instantiate so the compiler can generate them when compiling this cpp file
    template class AsyncLPRefiner<LocalKm1Policy>;
    template class AsyncLPRefiner<LocalCutPolicy>;

} // namespace mt_kahypar



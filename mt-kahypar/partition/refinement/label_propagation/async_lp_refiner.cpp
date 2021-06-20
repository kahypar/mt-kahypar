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
                                                     double) {
        ASSERT(_contraction_group_id != ds::invalidGroupID, "ContractionGroupID (Owner-ID) for locking is invalid.");
        ASSERT(!refinement_nodes.empty(), "AsyncLPRefiner will not work without given seed refinement nodes. Cannot be used "
                                          "solely for rebalancing or for global refinement!");
        ASSERT(std::all_of(refinement_nodes.begin(),refinement_nodes.end(),[&](const HypernodeID& hn) {return hypergraph.nodeIsEnabled(hn);})
               && "Not all given seed nodes are enabled!");

        _active_nodes.assign(refinement_nodes.begin(), refinement_nodes.end());

        _gain.reset();

        // Perform Label Propagation
        labelPropagation(hypergraph);

        // Update global part weight and sizes
        double imbalance;
        do {
            imbalance = metrics::imbalance(hypergraph, _context);
        } while (!best_metrics.update_imbalance_strong(imbalance));

        // Update metrics statistics
        Gain delta = _gain.delta();
        ASSERT(delta <= 0, "LP refiner worsen solution quality" << V(delta));
        best_metrics.fetch_add(delta, kahypar::Mode::direct_kway, _context.partition.objective);
        utils::Stats::instance().update_stat("lp_improvement", std::abs(delta));
        return delta < 0;

    }

    template <template <typename> class LocalGainPolicy>
    void AsyncLPRefiner<LocalGainPolicy>::labelPropagation(PartitionedHypergraph &hypergraph) {

        for (size_t i = 0; i < _context.refinement.label_propagation.maximum_iterations; ++i) {

            _next_active_nodes.clear();
            _visited_edges.clear();

            if ( !_active_nodes.empty() ) {
                labelPropagationRound(hypergraph, _next_active_nodes, _visited_edges);
            }

            _active_nodes = _next_active_nodes;

            auto has_node_duplicates = [&]() -> bool {
                std::sort(_next_active_nodes.begin(),_next_active_nodes.end());
                return std::adjacent_find(_next_active_nodes.begin(), _next_active_nodes.end()) != _next_active_nodes.end();
            };
            unused(has_node_duplicates);
            HEAVY_REFINEMENT_ASSERT(! has_node_duplicates());

            auto has_edge_duplicates = [&]() -> bool {
                std::sort(_visited_edges.begin(), _visited_edges.end());
                return std::adjacent_find(_visited_edges.begin(), _visited_edges.end()) != _visited_edges.end();
            };
            unused(has_edge_duplicates);
            HEAVY_REFINEMENT_ASSERT(! has_edge_duplicates());

            // Linear reset of anti-duplicator flags that were set in the last round
            for (const auto& hn : _next_active_nodes) {
                ASSERT(_next_active->isSet(hn), V(hn));
                bool reset = _next_active->compare_and_set_to_false(hn);
                unused(reset);
                ASSERT(reset);
            }
            for (const auto& he : _visited_edges) {
                ASSERT(_visited_he->isSet(he), V(he));
                bool reset = _visited_he->compare_and_set_to_false(he);
                unused(reset);
                ASSERT(reset);
            }

            if ( _active_nodes.empty() ) {
                break;
            }
        }
    }

    template <template <typename> class LocalGainPolicy>
    bool AsyncLPRefiner<LocalGainPolicy>::labelPropagationRound(PartitionedHypergraph &hypergraph,
                                                                NextActiveNodes &next_active_nodes,
                                                                VisitedEdges &visited_edges) {

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
            bool tried_to_move = false;
            // Lock node while moving to prevent concurrent uncontractions
            bool acquired = _lock_manager->tryToAcquireLock(hn, _contraction_group_id);
            if (acquired) {
              // Lock node in FM node tracker to make sure not to mess up PQs in any FM search
              bool locked_for_fm = _fm_node_tracker->tryAcquireNode(hn, _contraction_group_id);
              if (locked_for_fm) {
                tried_to_move = true;
                if ( ! moveVertex(hypergraph, hn, next_active_nodes, visited_edges, objective_delta) ) {
                  converged = false;
                }
                _fm_node_tracker->releaseNode(hn, _contraction_group_id);
              }
              _lock_manager->strongReleaseLock(hn, _contraction_group_id);
            }
            if (!tried_to_move) {
              retry_nodes.push_back(hn);
            }
        }

        // Retry acquiring lock and moving exactly once
        // todo mlaupichler: Number of retries is arbitrary, perhaps add CL option to set how often we retry here
        for ( size_t j = 0; j < retry_nodes.size(); ++j ) {
            const HypernodeID hn = retry_nodes[j];
            bool acquired = _lock_manager->tryToAcquireLock(hn, _contraction_group_id);
            if (acquired) {
              // Lock node in FM node tracker to make sure not to mess up PQs in any FM search
              bool locked_for_fm = _fm_node_tracker->tryAcquireNode(hn, _contraction_group_id);
              if (locked_for_fm) {
                if ( ! moveVertex(hypergraph, hn, next_active_nodes, visited_edges, objective_delta) ) {
                  converged = false;
                }
                _fm_node_tracker->releaseNode(hn, _contraction_group_id);
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



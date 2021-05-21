//
// Created by mlaupichler on 04.05.21.
//

#include <tbb/parallel_sort.h>
#include "asynch_lp_refiner.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

    template<template <typename> class LocalGainPolicy>
    bool AsynchLPRefiner<LocalGainPolicy>::refineImpl(PartitionedHypergraph &hypergraph,
                                                 const parallel::scalable_vector <mt_kahypar::HypernodeID> &refinement_nodes,
                                                 metrics::ThreadSafeMetrics &best_metrics, double) {
        ASSERT(_contraction_group_id != ds::invalidGroupID, "ContractionGroupID (Owner-ID) for locking is invalid.");
        ASSERT(!refinement_nodes.empty(), "AsynchLPRefiner will not work without given seed refinement nodes. Cannot be used "
                                          "solely for rebalancing or for global refinement!");
        ASSERT(std::all_of(refinement_nodes.begin(),refinement_nodes.end(),[&](const HypernodeID& hn) {return _lock_manager->isHeldBy(hn,_contraction_group_id);})
            && "Not all given seed nodes are locked by the contraction group id that this LP refinement is based on!");
        ASSERT(std::all_of(refinement_nodes.begin(),refinement_nodes.end(),[&](const HypernodeID& hn) {return hypergraph.nodeIsEnabled(hn);})
               && "Not all given seed nodes are enabled!");


        _seeds = refinement_nodes;
        _active_nodes.assign(refinement_nodes.begin(), refinement_nodes.end());

        _gain.reset();
        _next_active.reset();

        // Perform Label Propagation
        labelPropagation(hypergraph);

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

        HEAVY_REFINEMENT_ASSERT(hypergraph.checkTrackedPartitionInformation());
        HEAVY_REFINEMENT_ASSERT(current_metric + delta == metrics::objective(hypergraph, _context.partition.objective, false),
                                V(current_metric) << V(delta) <<
                                                  V(metrics::objective(hypergraph, _context.partition.objective, false)));

        best_metrics.fetch_add(delta, kahypar::Mode::direct_kway, _context.partition.objective);
        utils::Stats::instance().update_stat("lp_improvement", std::abs(delta));
        return delta < 0;

    }

    template <template <typename> class LocalGainPolicy>
    void AsynchLPRefiner<LocalGainPolicy>::labelPropagation(PartitionedHypergraph& hypergraph) {
        NextActiveNodes next_active_nodes;
        for (size_t i = 0; i < _context.refinement.label_propagation.maximum_iterations; ++i) {
            DBG << "Starting Label Propagation Round" << i;

            utils::Timer::instance().start_timer(
                    "lp_round_" + std::to_string(i), "Label Propagation Round " + std::to_string(i), true);

            if ( !_active_nodes.empty() ) {
                labelPropagationRound(hypergraph, next_active_nodes);
            }

            _active_nodes = next_active_nodes;
            next_active_nodes.clear();

            utils::Timer::instance().stop_timer("lp_round_" + std::to_string(i));

            if ( _active_nodes.empty() ) {
                break;
            }
        }

        // If the maximum iteration number was hit and thus active nodes remain, release their locks (except for the seeds)
        if (!_active_nodes.empty()) {
            auto seedIndices = getSortedIndicesOfSeedsInActiveNodes();
            size_t curSeed = 0;
            for ( size_t j = 0; j < _active_nodes.size(); ++j ) {
                ASSERT(_lock_manager->isHeldBy(_active_nodes[j], _contraction_group_id));
                // Skip seed indices
                if (curSeed < seedIndices.size() && j == seedIndices[curSeed]) {
                    ++curSeed;
                    continue;
                }
                const HypernodeID hn = _active_nodes[j];
                _lock_manager->strongReleaseLock(hn, _contraction_group_id);
            }
        }
    }

    template <template <typename> class LocalGainPolicy>
    bool AsynchLPRefiner<LocalGainPolicy>::labelPropagationRound(
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
        std::shuffle(_active_nodes.begin(),_active_nodes.end(), _rng);

        bool converged = true;
        for ( size_t j = 0; j < _active_nodes.size(); ++j ) {
            const HypernodeID hn = _active_nodes[j];
            if ( ! moveVertex(hypergraph, hn, next_active_nodes, objective_delta) ) {
                converged = false;
            }
        }

        // todo mlaupichler Releasing can probably be optimized. Releasing at the end of the round results in locks
        //  being held longer than they would optimally need to be held.

        // Release locks for nodes that were active in this round but will not be active in the next one. Exempt the
        // seed nodes of the refinement as their locks are never released by the refiner.
        // Do this only now at the end of the round in order to guarantee that during the round no locks will be released,
        // so queries whether a lock is already held by this LP hold until the end of the round.
        auto seedIndices = getSortedIndicesOfSeedsInActiveNodes();
        size_t curSeed = 0;
        for ( size_t j = 0; j < _active_nodes.size(); ++j ) {

            ASSERT(_lock_manager->isHeldBy(_active_nodes[j], _contraction_group_id));

            // Skip seed indices
            if (curSeed < seedIndices.size() && j == seedIndices[curSeed]) {
                ++curSeed;
                continue;
            }

            const HypernodeID hn = _active_nodes[j];
            if (!_next_active[hn]) {
                _lock_manager->strongReleaseLock(hn, _contraction_group_id);
            }
        }

        HEAVY_REFINEMENT_ASSERT(hypergraph.checkTrackedPartitionInformation());
        return converged;
    }

    template <template <typename> class LocalGainPolicy>
    parallel::scalable_vector<size_t> AsynchLPRefiner<LocalGainPolicy>::getSortedIndicesOfSeedsInActiveNodes() const {
        parallel::scalable_vector<size_t> indices;

        for ( size_t j = 0; j < _active_nodes.size(); ++j ) {
          // REVIEW should be fine in most cases, but looks dangerous.
          // in case _seeds.size() is large, e.g., > 40 or smth, you can sort the seeds and do binary search (std::lower_bound) whether _active_nodes[j] is a seed
          // also, you can do this directly in the loop above. no need to allocate a vector
            for (HypernodeID seed : _seeds) {
                if (seed == _active_nodes[j]) indices.push_back(j);
            }
        }
        ASSERT(indices.size() <= _seeds.size());
        HEAVY_REFINEMENT_ASSERT(std::is_sorted(indices.begin(),indices.end()) && "Seed indices are not sorted!");
        return indices;
    }

    template <template <typename> class LocalGainPolicy>
    void AsynchLPRefiner<LocalGainPolicy>::resetForGroup(ds::ContractionGroupID groupID) {
        _contraction_group_id = groupID;
    }

    // explicitly instantiate so the compiler can generate them when compiling this cpp file
    template class AsynchLPRefiner<LocalKm1Policy>;
    template class AsynchLPRefiner<LocalCutPolicy>;

} // namespace mt_kahypar



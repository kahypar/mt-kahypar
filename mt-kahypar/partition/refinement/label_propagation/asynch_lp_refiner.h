//
// Created by mlaupichler on 04.05.21.
//

#ifndef KAHYPAR_ASYNCH_LP_REFINER_H
#define KAHYPAR_ASYNCH_LP_REFINER_H

#include <mt-kahypar/partition/refinement/i_refiner.h>
#include <mt-kahypar/partition/refinement/policies/gain_policy.h>

namespace mt_kahypar {

    /// Label Propagation Refiner to be used in asynchronous uncoarsening for localized refinement. Uses global locking
    /// datastructure to allow multiple concurrent localized label propagations as well as concurrent uncontractions.
    /// Can not be used for global refinement or (only) rebalancing as it always requires seed nodes for refining!
    template <template <typename> class GainPolicy> class AsynchLPRefiner : public IRefiner {
    private:
        using GainCalculator = GainPolicy<PartitionedHypergraph>;
        using ActiveNodes = parallel::scalable_vector<HypernodeID>;
        using NextActiveNodes = ds::StreamingVector<HypernodeID>;

        static constexpr bool debug = false;
        static constexpr bool enable_heavy_assert = false;

    public:
        explicit AsynchLPRefiner(Hypergraph &hypergraph, const Context &context, const TaskGroupID task_group_id,
                                 ds::IGroupLockManager *lockManager, ds::ContractionGroupID contraction_group_id) :
        _context(context),
        _task_group_id(task_group_id),
        _gain(context),
        _active_nodes(),
        _next_active(hypergraph.initialNumNodes()),
        _visited_he(hypergraph.initialNumEdges()),
        _lock_manager(lockManager),
        _contraction_group_id(contraction_group_id),
        _seeds(),
        _num_acquired_locks(0),
        _num_released_locks(0),
        _num_nodes(hypergraph.initialNumNodes()) { }

        AsynchLPRefiner(const AsynchLPRefiner&) = delete;
        AsynchLPRefiner(AsynchLPRefiner&&) = delete;

        AsynchLPRefiner & operator= (const AsynchLPRefiner &) = delete;
        AsynchLPRefiner & operator= (AsynchLPRefiner &&) = delete;

    private:

        bool refineImpl(PartitionedHypergraph& hypergraph,
                        const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                        kahypar::Metrics& best_metrics,
                        double) final ;

        void labelPropagation(PartitionedHypergraph& hypergraph);

        bool labelPropagationRound(PartitionedHypergraph& hypergraph, NextActiveNodes& next_active_nodes);

        template<typename F>
        bool moveVertex(PartitionedHypergraph& hypergraph,
                        const HypernodeID hn,
                        NextActiveNodes& next_active_nodes,
                        const F& objective_delta) {
            bool is_moved = false;
            ASSERT(hn != kInvalidHypernode);
            if ( hypergraph.isBorderNode(hn) ) {
                ASSERT(hypergraph.nodeIsEnabled(hn));

                Move best_move = _gain.computeMaxGainMove(hypergraph, hn);
                // We perform a move if it either improves the solution quality or, in case of a
                // zero gain move, the balance of the solution.
                const bool positive_gain = best_move.gain < 0;
                const bool zero_gain_move = (_context.refinement.label_propagation.rebalancing &&
                                             best_move.gain == 0 &&
                                             hypergraph.partWeight(best_move.from) - 1 >
                                             hypergraph.partWeight(best_move.to) + 1 &&
                                             hypergraph.partWeight(best_move.to) <
                                             _context.partition.perfect_balance_part_weights[best_move.to]);
                const bool perform_move = positive_gain || zero_gain_move;
                if (best_move.from != best_move.to && perform_move) {
                    PartitionID from = best_move.from;
                    PartitionID to = best_move.to;

                    Gain delta_before = _gain.localDelta();
                    bool changed_part = changeNodePart(hypergraph, hn, from, to, objective_delta);
                    if (changed_part) {
                        // In case the move to block 'to' was successful, we verify that the "real" gain
                        // of the move is either equal to our computed gain or if not, still improves
                        // the solution quality.
                        Gain move_delta = _gain.localDelta() - delta_before;
                        bool accept_move = (move_delta == best_move.gain || move_delta <= 0);
                        if (accept_move) {
                            DBG << "Move hypernode" << hn << "from block" << from << "to block" << to
                                << "with gain" << best_move.gain << "( Real Gain: " << move_delta << ")";

                            // Set all neighbors of the vertex to active
                            for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
                                if ( hypergraph.edgeSize(he) <=
                                     ID(_context.refinement.label_propagation.hyperedge_size_activation_threshold) ) {
                                    if ( !_visited_he[he] ) {
                                        for (const HypernodeID& pin : hypergraph.pins(he)) {
                                            if (!_lock_manager->isLocked(pin)) {
                                                bool acquired = _lock_manager->tryToAcquireLock(pin, _contraction_group_id);
                                                // If lock could be acquired this has to be the first time that this hypernode has been traversed in this LP round
                                                // => Assert that it is inserted into the next active nodes
                                                if (acquired) {
                                                    ++_num_acquired_locks;
                                                    bool set_next_active = _next_active.compare_and_set_to_true(pin);
                                                    ASSERT(set_next_active);
                                                    next_active_nodes.stream(pin);
                                                }
                                            } else if (_lock_manager->isHeldBy(pin,_contraction_group_id) && _next_active.compare_and_set_to_true(pin) ) {
                                                next_active_nodes.stream(pin);
                                            }
                                        }
                                        _visited_he.set(he, true);
                                    }
                                }
                            }
                            if ( _next_active.compare_and_set_to_true(hn) ) {
                                ASSERT(_lock_manager->isHeldBy(hn,_contraction_group_id));
                                next_active_nodes.stream(hn);
                            }
                            is_moved = true;
                        } else {
                            DBG << "Revert move of hypernode" << hn << "from block" << from << "to block" << to
                                << "( Expected Gain:" << best_move.gain << ", Real Gain:" << move_delta << ")";
                            // In case, the real gain is not equal to the computed gain and
                            // worsens the solution quality we revert the move.
                            ASSERT(hypergraph.partID(hn) == to);
                            changeNodePart(hypergraph, hn, to, from, objective_delta);
                        }
                    }
                }
            }

            return is_moved;
        }

        // NOOP as the asynch refiner is never supposed to be used for global refinement.
        // If initialize is called and then refine is called without giving refinement nodes as parameters,
        // the assertion in refineImpl will fail.
        void initializeImpl(PartitionedHypergraph&) final {};

        template<typename F>
        bool changeNodePart(PartitionedHypergraph& phg,
                            const HypernodeID hn,
                            const PartitionID from,
                            const PartitionID to,
                            const F& objective_delta) {
            bool success = false;
            if ( _context.partition.paradigm == Paradigm::nlevel && phg.isGainCacheInitialized()) {
                success = phg.changeNodePartWithGainCacheUpdate(hn, from, to,
                                                                _context.partition.max_part_weights[to], [] { }, objective_delta);
            } else {
                success = phg.changeNodePart(hn, from, to,
                                             _context.partition.max_part_weights[to], []{}, objective_delta);
            }
            return success;
        }

        parallel::scalable_vector<size_t> getSortedIndicesOfSeedsInActiveNodes() const;

        const Context& _context;
        const TaskGroupID _task_group_id;
        GainCalculator _gain;
        ActiveNodes _active_nodes;
        ds::ThreadSafeFastResetFlagArray<> _next_active;
        kahypar::ds::FastResetFlagArray<> _visited_he;

        // todo mlaupichler The AsynchLPRefiner should not be bound so tightly to the idea of ContractionGroups. Use a template argument for the OwnerID (also in the _lock_manager pointer) in order to allow any OwnerID for the locking.
        // ! The ID of the contraction group that the seed refinement nodes are from. Used as an identifier for locking
        // nodes so uncontracting and refinement of a contraction group use the same locks.
        ds::ContractionGroupID _contraction_group_id;

        // ! A pointer to the LockManager that is used by uncontraction and refinement operations during the entire
        // uncoarsening.
        ds::IGroupLockManager* _lock_manager;

        // ! A reference to the seed nodes for the current refinement. Locks for seed nodes are never
        parallel::scalable_vector<HypernodeID> _seeds;

        // todo mlaupichler remove debug: these two counters
        size_t _num_acquired_locks;
        size_t _num_released_locks;
        size_t _num_nodes;
    };

    using AsynchLPKm1Refiner = AsynchLPRefiner<Km1Policy>;
    using AsynchLPCutRefiner = AsynchLPRefiner<CutPolicy>;
}  // namespace kahypar



#endif //KAHYPAR_ASYNCH_LP_REFINER_H

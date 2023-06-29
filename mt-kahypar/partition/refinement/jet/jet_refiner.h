/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Nikolai Maas <nikolai.maas@kit.edu>
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#pragma once

#include "kahypar/datastructure/fast_reset_flag_array.h"

#include "mt-kahypar/datastructures/thread_safe_fast_reset_flag_array.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/utils/cast.h"


namespace mt_kahypar {
template <typename TypeTraits, typename GainTypes>
class JetRefiner final : public IRefiner {
 private:
  using Hypergraph = typename TypeTraits::Hypergraph;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
  using GainCache = typename GainTypes::GainCache;
  using GainCalculator = typename GainTypes::GainComputation;
  using RatingMap = typename GainCalculator::RatingMap;
  using AttributedGains = typename GainTypes::AttributedGains;
  using ActiveNodes = parallel::scalable_vector<HypernodeID>;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

 public:
  explicit JetRefiner(const HypernodeID num_hypernodes,
                      const HyperedgeID num_hyperedges,
                      const Context& context,
                      GainCache& gain_cache,
                      IRefiner& rebalancer);

  explicit JetRefiner(const HypernodeID num_hypernodes,
                      const HyperedgeID num_hyperedges,
                      const Context& context,
                      gain_cache_t gain_cache,
                      IRefiner& rebalancer) :
    JetRefiner(num_hypernodes, num_hyperedges, context,
               GainCachePtr::cast<GainCache>(gain_cache), rebalancer) {}

  JetRefiner(const JetRefiner&) = delete;
  JetRefiner(JetRefiner&&) = delete;

  JetRefiner & operator= (const JetRefiner &) = delete;
  JetRefiner & operator= (JetRefiner &&) = delete;

 private:
  bool refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                  const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                  Metrics& best_metrics,
                  double time_limit) final;

  void labelPropagationRound(PartitionedHypergraph& hypergraph);

  template<typename F>
  bool moveVertexGreedily(PartitionedHypergraph& hypergraph,
                          const HypernodeID hn,
                          const F& objective_delta) {
    bool is_moved = false;
    ASSERT(hn != kInvalidHypernode);
    if ( hypergraph.isBorderNode(hn) ) {
      ASSERT(hypergraph.nodeIsEnabled(hn));

      Move best_move;
      if (_context.refinement.jet.unconstrained_upper_bound >= 1.0) {
        best_move = _gain.computeMaxGainMove(hypergraph, hn, false, false, false,
                                             _context.refinement.jet.unconstrained_upper_bound);
      } else {
        best_move = _gain.computeMaxGainMove(hypergraph, hn, false, false, true);
      }
      const bool positive_gain = best_move.gain < 0;
      if (positive_gain && best_move.from != best_move.to) {
        PartitionID from = best_move.from;
        PartitionID to = best_move.to;

        Gain delta_before = _gain.localDelta();
        if (changeNodePart(hypergraph, hn, from, to, objective_delta)) {
          is_moved = true;

          // In case the move to block 'to' was successful, we verify that the "real" gain
          // of the move is either equal to our computed gain or if not, still improves
          // the solution quality.
          Gain move_delta = _gain.localDelta() - delta_before;
          bool accept_move = (move_delta == best_move.gain || move_delta <= 0);
          if (!accept_move) {
            ASSERT(hypergraph.partID(hn) == to);
            changeNodePart(hypergraph, hn, to, from, objective_delta);
          }
        }
      }
    }

    return is_moved;
  }

  void initializeActiveNodes(PartitionedHypergraph& hypergraph,
                             const parallel::scalable_vector<HypernodeID>& refinement_nodes);

  void initializeImpl(mt_kahypar_partitioned_hypergraph_t& phg) final;

  template<bool precomputed>
  void computeActiveNodesFromGraph(const PartitionedHypergraph& hypergraph);

  void computeActiveNodesFromVector(const PartitionedHypergraph& hypergraph,
                                    const parallel::scalable_vector<HypernodeID>& refinement_nodes);

  template <bool precomputed>
  void computeActiveNodesFromPreviousRound(const PartitionedHypergraph& hypergraph);

  // ! Applied during computation of active nodes. If precomputed, applies the first JET filter
  template<bool precomputed, typename F>
  void processNode(const PartitionedHypergraph& hypergraph, const HypernodeID hn, F add_node_fn, const bool top_level);

  void recomputePenalties(const PartitionedHypergraph& hypergraph, bool did_rebalance);

  void storeCurrentPartition(const PartitionedHypergraph& hypergraph, parallel::scalable_vector<PartitionID>& parts);

  void rollbackToBestPartition(PartitionedHypergraph& hypergraph);

  void rebalance(PartitionedHypergraph& hypergraph, Metrics& current_metrics, double time_limit);

  template<typename F>
  bool changeNodePart(PartitionedHypergraph& phg,
                      const HypernodeID hn,
                      const PartitionID from,
                      const PartitionID to,
                      const F& objective_delta,
                      bool force_move = false) {
    HypernodeWeight max_weight;
    if (_context.refinement.jet.unconstrained_upper_bound >= 1.0 && !force_move) {
      max_weight = _context.refinement.jet.unconstrained_upper_bound * _context.partition.max_part_weights[to];
    } else {
      max_weight = std::numeric_limits<HypernodeWeight>::max();
    }
    if ( _context.forceGainCacheUpdates() && _gain_cache.isInitialized() ) {
      return phg.changeNodePart(_gain_cache, hn, from, to, max_weight, []{}, objective_delta);
    } else {
      return phg.changeNodePart(hn, from, to, max_weight, []{}, objective_delta);
    }
  }

  void resizeDataStructuresForCurrentK() {
    // If the number of blocks changes, we resize data structures
    // (can happen during deep multilevel partitioning)
    if ( _current_k != _context.partition.k ) {
      _current_k = _context.partition.k;
      _gain.changeNumberOfBlocks(_current_k);
    }
  }

  const Context& _context;
  GainCache& _gain_cache;
  bool _precomputed;
  PartitionID _current_k;
  HypernodeID _top_level_num_nodes;
  bool _current_partition_is_best;
  parallel::scalable_vector<PartitionID> _best_partition;
  parallel::scalable_vector<PartitionID> _current_partition;
  GainCalculator _gain;
  ActiveNodes _active_nodes;
  parallel::scalable_vector<std::pair<Gain, PartitionID>> _gains_and_target;
  ds::ThreadSafeFastResetFlagArray<> _next_active;
  kahypar::ds::FastResetFlagArray<> _visited_he;
  IRefiner& _rebalancer;
};

}  // namespace kahypar

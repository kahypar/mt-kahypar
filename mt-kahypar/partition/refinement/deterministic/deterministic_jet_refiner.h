/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/datastructures/buffered_vector.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/i_rebalancer.h"

#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/utils/reproducible_random.h"

namespace mt_kahypar {

template<typename GraphAndGainTypes>
class DeterministicJetRefiner final : public IRefiner {

  using PartitionedHypergraph = typename GraphAndGainTypes::PartitionedHypergraph;
  using GainComputation = typename GraphAndGainTypes::GainComputation;
  using AttributedGains = typename GraphAndGainTypes::AttributedGains;
  using GainCache = typename GraphAndGainTypes::GainCache;
  using ActiveNodes = typename parallel::scalable_vector<HypernodeID>;
  using RatingMap = typename GainComputation::RatingMap;

public:

  explicit DeterministicJetRefiner(const HypernodeID num_hypernodes,
    const HyperedgeID num_hyperedges,
    const Context& context,
    gain_cache_t gain_cache,
    IRebalancer& rebalancer) :
    DeterministicJetRefiner(num_hypernodes, num_hyperedges, context,
      GainCachePtr::cast<GainCache>(gain_cache), rebalancer) {}

  explicit DeterministicJetRefiner(const HypernodeID num_hypernodes,
    const HyperedgeID,
    const Context& context,
    GainCache& gain_cache,
    IRebalancer& rebalancer) :
    _context(context),
    _current_k(context.partition.k),
    _top_level_num_nodes(num_hypernodes),
    _current_partition_is_best(true),
    _active_nodes(),
    _moves(),
    _best_partition(num_hypernodes, kInvalidPartition),
    _current_partition(num_hypernodes, kInvalidPartition),
    _gain_computation(context, true /* disable_randomization */),
    _gains_and_target(num_hypernodes),
    _locks(num_hypernodes),
    _rebalancer(rebalancer),
    _gain_cache(gain_cache) {}

private:
  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  bool refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
    const vec<HypernodeID>&,
    Metrics& best_metrics, double) final;

  void initializeImpl(mt_kahypar_partitioned_hypergraph_t& phg);

  void computeActiveNodesFromGraph(const PartitionedHypergraph& hypergraph, bool first_round);

  Gain performMoveWithAttributedGain(PartitionedHypergraph& phg, const HypernodeID hn);

  void storeCurrentPartition(const PartitionedHypergraph& hypergraph, parallel::scalable_vector<PartitionID>& parts);

  void rollbackToBestPartition(PartitionedHypergraph& hypergraph);

  template<typename F>
  void changeNodePart(PartitionedHypergraph& phg,
    const HypernodeID hn,
    const PartitionID from,
    const PartitionID to,
    const F& objective_delta) {
    bool success = false;
    constexpr HypernodeWeight inf_weight = std::numeric_limits<HypernodeWeight>::max();
    if (_gain_cache.isInitialized()) {
      success = phg.changeNodePart(_gain_cache, hn, from, to, inf_weight, [] {}, objective_delta);
    } else {
      success = phg.changeNodePart(hn, from, to, inf_weight, [] {}, objective_delta);
    }
    ASSERT(success);
    unused(success);
  }

  void recomputePenalties(const PartitionedHypergraph& hypergraph, bool did_rebalance);

  bool arePotentialMovesToOtherParts(const PartitionedHypergraph& hypergraph, const parallel::scalable_vector<HypernodeID>& moves);

  bool noInvalidPartitions(const PartitionedHypergraph& phg, const parallel::scalable_vector<PartitionID>& parts);

  void resizeDataStructuresForCurrentK() {
    // If the number of blocks changes, we resize data structures
    // (can happen during deep multilevel partitioning)
    if (_current_k != _context.partition.k) {
      _current_k = _context.partition.k;
      _gain_computation.changeNumberOfBlocks(_current_k);
      if (_gain_cache.isInitialized()) {
        _gain_cache.changeNumberOfBlocks(_current_k);
      }
    }
  }

  const Context& _context;
  PartitionID _current_k;
  HypernodeID _top_level_num_nodes;
  bool _current_partition_is_best;
  ActiveNodes _active_nodes;
  parallel::scalable_vector<HypernodeID> _moves;
  parallel::scalable_vector<PartitionID> _best_partition;
  parallel::scalable_vector<PartitionID> _current_partition;
  GainComputation _gain_computation;
  parallel::scalable_vector<std::pair<Gain, PartitionID>> _gains_and_target;
  kahypar::ds::FastResetFlagArray<> _locks;
  IRebalancer& _rebalancer;
  GainCache& _gain_cache;

};

}
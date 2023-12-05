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
    GainCache&,
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
    tmp_active_nodes(),
    _part_before_round(num_hypernodes),
    _afterburner_gain(num_hypernodes),
    _afterburner_buffer(_current_k, 0),
    _hyperedge_buffer() {}

private:
  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  bool refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
    const vec<HypernodeID>&,
    Metrics& best_metrics, double) final;

  void initializeImpl(mt_kahypar_partitioned_hypergraph_t& phg);

  void computeActiveNodesFromGraph(const PartitionedHypergraph& hypergraph);

  Gain performMoveWithAttributedGain(PartitionedHypergraph& phg, const HypernodeID hn);

  void storeCurrentPartition(const PartitionedHypergraph& hypergraph, parallel::scalable_vector<PartitionID>& parts);

  void rollbackToBestPartition(PartitionedHypergraph& hypergraph);

  void changeNodePart(PartitionedHypergraph& phg,
    const HypernodeID hn,
    const PartitionID from,
    const PartitionID to) {
    constexpr HypernodeWeight inf_weight = std::numeric_limits<HypernodeWeight>::max();
    const bool success = phg.changeNodePartNoSync(hn, from, to, inf_weight);
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
    }
  }

  HyperedgeWeight calculateGainDelta(const PartitionedHypergraph& phg) const {
    tbb::enumerable_thread_specific<HyperedgeWeight> gain_delta(0);
    phg.doParallelForAllNodes([&](const HypernodeID hn) {
      const PartitionID from = _part_before_round[hn];
      const PartitionID to = phg.partID(hn);
      if (from != to) {
        for (const HyperedgeID& he : phg.incidentEdges(hn)) {
          HypernodeID pin_count_in_from_part_after = 0;
          HypernodeID pin_count_in_to_part_after = 1;
          for (const HypernodeID& pin : phg.pins(he)) {
            if (pin != hn) {
              const PartitionID part = pin < hn ? phg.partID(pin) : _part_before_round[pin];
              if (part == from) {
                pin_count_in_from_part_after++;
              } else if (part == to) {
                pin_count_in_to_part_after++;
              }
            }
          }
          SynchronizedEdgeUpdate sync_update;
          sync_update.he = he;
          sync_update.edge_weight = phg.edgeWeight(he);
          sync_update.edge_size = phg.edgeSize(he);
          sync_update.pin_count_in_from_part_after = pin_count_in_from_part_after;
          sync_update.pin_count_in_to_part_after = pin_count_in_to_part_after;
          gain_delta.local() += AttributedGains::gain(sync_update);
        }
      }
    });
    return gain_delta.combine(std::plus<>());
  }

  void hypergraphAfterburner(const PartitionedHypergraph& phg) {
    tbb::parallel_for(0UL, _afterburner_gain.size(), [&](const size_t i) {
      _afterburner_gain[i].store(0);
    });
    phg.doParallelForAllEdges([&](const HyperedgeID& he) {
      auto& edgeBuffer = _hyperedge_buffer.local();
      auto& afterburnerBuffer = _afterburner_buffer.local();
      const HypernodeID edgeSize = phg.edgeSize(he);
      if (edgeSize > edgeBuffer.size()) {
        edgeBuffer.resize(edgeSize);
      }
      // materialize Hyperedge
      size_t index = 0;
      for (const auto pin : phg.pins(he)) {
        edgeBuffer[index] = pin;
        ++index;
      }
      // sort by afterburner order
      std::sort(edgeBuffer.begin(), edgeBuffer.begin() + edgeSize, [&](const HypernodeID& a, const HypernodeID& b) {
        auto [gain_a, to_a] = _gains_and_target[a];
        auto [gain_b, to_b] = _gains_and_target[b];
        return (gain_a < gain_b || (gain_a == gain_b && a < b));
      });

      // initial pin-counts
      for (auto& pinCount : afterburnerBuffer) {
        pinCount = 0;
      }
      for (size_t i = 0; i < edgeSize; ++i) {
        const HypernodeID pin = edgeBuffer[i];
        afterburnerBuffer[phg.partID(pin)]++;
      }
      // update pin-counts for each pin
      for (size_t i = 0; i < edgeSize; ++i) {
        const HypernodeID pin = edgeBuffer[i];
        const PartitionID from = phg.partID(pin);
        const auto [gain, to] = _gains_and_target[pin];
        afterburnerBuffer[from]--;
        afterburnerBuffer[to]++;
        SynchronizedEdgeUpdate sync_update;
        sync_update.he = he;
        sync_update.edge_weight = phg.edgeWeight(he);
        sync_update.edge_size = phg.edgeSize(he);
        sync_update.pin_count_in_from_part_after = afterburnerBuffer[from];
        sync_update.pin_count_in_to_part_after = afterburnerBuffer[to];
        const Gain attributedGain = AttributedGains::gain(sync_update);
        if (gain != 0) {
          _afterburner_gain[pin] += attributedGain;
        }
      }
    });
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
  ds::StreamingVector<HypernodeID> tmp_active_nodes;
  parallel::scalable_vector<PartitionID> _part_before_round;

  // hypergraph afterburner
  parallel::scalable_vector<std::atomic<Gain>> _afterburner_gain;
  tbb::enumerable_thread_specific<std::vector<size_t>> _afterburner_buffer;
  tbb::enumerable_thread_specific<std::vector<HypernodeID>> _hyperedge_buffer;
};

}

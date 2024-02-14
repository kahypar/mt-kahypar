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
    const HyperedgeID num_hyperedges,
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
    _hyperedge_buffer(),
    _edge_flag(num_hyperedges),
    _current_edge_flag(1) {}

private:
  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  bool refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
    const vec<HypernodeID>&,
    Metrics& best_metrics, double) final;

  void initializeImpl(mt_kahypar_partitioned_hypergraph_t& phg);

  void computeActiveNodesFromGraph(const PartitionedHypergraph& hypergraph);

  void storeCurrentPartition(const PartitionedHypergraph& hypergraph, parallel::scalable_vector<PartitionID>& parts);

  template <bool isGraph>
  void rollbackToBestPartition(PartitionedHypergraph& phg) {
    // This function is passed as lambda to the changeNodePart function and used
    // to calculate the "real" delta of a move (in terms of the used objective function).
    auto objective_delta = [&](const SynchronizedEdgeUpdate& sync_update) {
      _gain_computation.computeDeltaForHyperedge(sync_update);
    };

    auto reset_node = [&](const HypernodeID hn) {
      const PartitionID part_id = phg.partID(hn);
      if (part_id != _best_partition[hn]) {
        ASSERT(_best_partition[hn] != kInvalidPartition);
        ASSERT(_best_partition[hn] >= 0 && _best_partition[hn] < _current_k);
        changeNodePart<isGraph>(phg, hn, part_id, _best_partition[hn], objective_delta);
      }
    };
    phg.doParallelForAllNodes(reset_node);
    _current_partition_is_best = true;
  }

  template< bool isGraph, typename F>
  void changeNodePart(PartitionedHypergraph& phg,
    const HypernodeID hn,
    const PartitionID from,
    const PartitionID to,
    const F& objective_delta) {
    constexpr HypernodeWeight inf_weight = std::numeric_limits<HypernodeWeight>::max();
    const bool success = isGraph ? phg.changeNodePartNoSync(hn, from, to, inf_weight) : phg.changeNodePart(hn, from, to, inf_weight, [] {}, objective_delta);
    ASSERT(success);
    unused(success);
  }

  template<bool isGraph>
  Gain performMoveWithAttributedGain(
    PartitionedHypergraph& phg, const HypernodeID hn) {
    const auto from = phg.partID(hn);
    const auto [gain, to] = _gains_and_target[hn];
    Gain attributed_gain = 0;
    auto objective_delta = [&](const SynchronizedEdgeUpdate& sync_update) {
      attributed_gain -= AttributedGains::gain(sync_update);
    };
    ASSERT(to >= 0 && to < _current_k);
    changeNodePart<isGraph>(phg, hn, from, to, objective_delta);
    return attributed_gain;
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

    auto hardcoded_afterburn = [&](const HyperedgeID& he) {
      HypernodeID a = kInvalidHypernode;
      HypernodeID b = kInvalidHypernode;
      Gain minGain = std::numeric_limits<Gain>::max();
      for (const auto pin : phg.pins(he)) {
        const Gain& gain = _gains_and_target[pin].first;
        if (a == kInvalidHypernode || gain < minGain || (minGain == gain && pin < a)) {
          b = a;
          a = pin;
          minGain = gain;
        } else {
          b = pin;
        }
      }
      const auto& [gain_a, to_a] = _gains_and_target[a];
      const auto& [gain_b, to_b] = _gains_and_target[b];
      const PartitionID from_a = phg.partID(a);
      const PartitionID from_b = phg.partID(b);
      const HyperedgeWeight weight = phg.edgeWeight(he);
      // moving a
      if (from_a != to_a) {
        if (from_a == from_b) {
          _afterburner_gain[a] += weight;
        } else if (to_a == from_b) {
          _afterburner_gain[a] -= weight;
        }
      }
      //moving b after a
      if (from_b != to_b) {
        if (from_b == to_a) {
          _afterburner_gain[b] += weight;
        } else if (to_b == to_a) {
          _afterburner_gain[b] -= weight;
        }
      }
    };

    auto afterburn_edge = [&](const HyperedgeID& he) {
      const HypernodeID edgeSize = phg.edgeSize(he);
      if (_context.refinement.deterministic_refinement.jet.afterburner_hardcode_graph_edges && edgeSize == 2) return hardcoded_afterburn(he);
      auto& edgeBuffer = _hyperedge_buffer.local();
      auto& afterburnerBuffer = _afterburner_buffer.local();
      if (edgeSize > edgeBuffer.size()) {
        edgeBuffer.resize(edgeSize);
      }
      // initial pin-counts
      for (auto& pinCount : afterburnerBuffer) {
        pinCount = 0;
      }

      // materialize Hyperedge
      size_t index = 0;
      for (const auto pin : phg.pins(he)) {
        const auto part = phg.partID(pin);
        if (!_context.refinement.deterministic_refinement.jet.afterburner_skip_unmoved_pins || part != _gains_and_target[pin].second) {
          edgeBuffer[index] = pin;
          ++index;
        } else {
          afterburnerBuffer[part]++;
        }
      }
      if (_context.refinement.deterministic_refinement.jet.afterburner_sorting_nets && index < 4) {
        if (index == 0) {
          return;
        } else if (index == 1) {
        } else if (index == 2) {
          auto& a = edgeBuffer[0];
          auto& b = edgeBuffer[1];
          auto& [gain_a, to_a] = _gains_and_target[a];
          auto& [gain_b, to_b] = _gains_and_target[b];
          if (!(gain_a < gain_b || (gain_a == gain_b && a < b))) std::swap(a, b);
        } else if (index == 3) {
          auto& a = edgeBuffer[0];
          auto& b = edgeBuffer[1];
          auto& c = edgeBuffer[2];
          if (!(_gains_and_target[a].first < _gains_and_target[c].first || (_gains_and_target[a].first == _gains_and_target[c].first && a < c)))
            std::swap(a, c);
          if (!(_gains_and_target[a].first < _gains_and_target[b].first || (_gains_and_target[a].first == _gains_and_target[b].first && a < b)))
            std::swap(a, b);
          if (!(_gains_and_target[b].first < _gains_and_target[c].first || (_gains_and_target[b].first == _gains_and_target[c].first && b < c)))
            std::swap(b, c);
        }
      } else {
        // sort by afterburner order
        std::sort(edgeBuffer.begin(), edgeBuffer.begin() + index, [&](const HypernodeID& a, const HypernodeID& b) {
          auto& [gain_a, to_a] = _gains_and_target[a];
          auto& [gain_b, to_b] = _gains_and_target[b];
          return (gain_a < gain_b || (gain_a == gain_b && a < b));
        });
      }
      for (size_t i = 0; i < index; ++i) {
        const HypernodeID pin = edgeBuffer[i];
        afterburnerBuffer[phg.partID(pin)]++;     // NOTE can this not be in the same loop that insets into edgeBuffer (L250), i.e., move L256 out of the else branch to always execute
      }
      // update pin-counts for each pin
      for (size_t i = 0; i < index; ++i) {
        const HypernodeID pin = edgeBuffer[i];
        const PartitionID from = phg.partID(pin);
        const auto [gain, to] = _gains_and_target[pin];
        afterburnerBuffer[from]--;
        afterburnerBuffer[to]++;
        SynchronizedEdgeUpdate sync_update;
        sync_update.he = he;                    // NOTE .he and .edge_weight inits can be moved outside the loop?
        sync_update.edge_weight = phg.edgeWeight(he);
        sync_update.edge_size = phg.edgeSize(he);
        sync_update.pin_count_in_from_part_after = afterburnerBuffer[from];
        sync_update.pin_count_in_to_part_after = afterburnerBuffer[to];
        const Gain attributedGain = AttributedGains::gain(sync_update);
        if (!_context.refinement.deterministic_refinement.jet.afterburner_skip_zero || attributedGain != 0) {
          _afterburner_gain[pin] += attributedGain;   // NOTE specify memory_order_relaxed (might be faster). Also beware that (in our twisted definition) positive attributed gain constitutes an increase in the objective function
        }
      }
    };

    if (_context.refinement.deterministic_refinement.jet.afterburner_incident_edges) {
      tbb::parallel_for(0UL, _active_nodes.size(), [&](const size_t& i) {
        const HypernodeID hn = _active_nodes[i];
        for (const HyperedgeID& he : phg.incidentEdges(hn)) {
          size_t flag = _edge_flag[he].load();    // NOTE specify memory_order_acquire in compare_exchange and memory_order_relaxed in load
          if (flag == _current_edge_flag || !_edge_flag[he].compare_exchange_strong(flag, _current_edge_flag)) continue;
          afterburn_edge(he);
        }
      });
      _current_edge_flag++;
    } else {
      phg.doParallelForAllEdges([&](const HyperedgeID& he) {
        afterburn_edge(he);
      });
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
  ds::StreamingVector<HypernodeID> tmp_active_nodes;
  parallel::scalable_vector<PartitionID> _part_before_round;

  // hypergraph afterburner
  parallel::scalable_vector<std::atomic<Gain>> _afterburner_gain;
  tbb::enumerable_thread_specific<std::vector<size_t>> _afterburner_buffer;
  tbb::enumerable_thread_specific<std::vector<HypernodeID>> _hyperedge_buffer;  // NOTE merge the two ETS into one with a struct --> only one lookup
  // incident edges in hypergraph afterburner
  parallel::scalable_vector<std::atomic<size_t>> _edge_flag;  // NOTE 16 bit should be enough
  size_t _current_edge_flag;
  double _negative_gain_factor;
};

}

/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
******************************************************************************/

#pragma once

#include "tbb/task.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/initial_partitioning/flat/initial_partitioning_data_container.h"

namespace mt_kahypar {
template<typename TypeTraits,
         typename GainPolicy,
         typename PQSelectionPolicy>
class GreedyInitialPartitionerT : public tbb::task {
  using HyperGraph = typename TypeTraits::HyperGraph;
  using InitialPartitioningDataContainer = InitialPartitioningDataContainerT<TypeTraits>;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;
  static PartitionID kInvalidPartition;
  static HypernodeID kInvalidHypernode;
  static Gain kInvalidGain;

 public:
  GreedyInitialPartitionerT(const InitialPartitioningAlgorithm algorithm,
                            InitialPartitioningDataContainer& ip_data,
                            const Context& context) :
    _algorithm(algorithm),
    _ip_data(ip_data),
    _context(context) { }

  tbb::task* execute() override {
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    HyperGraph& hg = _ip_data.local_hypergraph();
    KWayPriorityQueue& kway_pq = _ip_data.local_kway_priority_queue();
    kahypar::ds::FastResetFlagArray<>& hyperedges_in_queue =
      _ip_data.local_hyperedge_fast_reset_flag_array();

    // Insert start vertices into its corresponding PQs
    parallel::scalable_vector<HypernodeID> start_nodes =
      PseudoPeripheralStartNodes<TypeTraits>::computeStartNodes(_ip_data, _context);
    ASSERT(static_cast<size_t>(_context.partition.k) == start_nodes.size());
    kway_pq.clear();
    for ( PartitionID block = 0; block < _context.partition.k; ++block ) {
      insertVertexIntoPQ(hg, kway_pq, start_nodes[block], block);
    }

    _ip_data.reset_unassigned_hypernodes();
    hyperedges_in_queue.reset();
    PartitionID to = kInvalidPartition;
    bool use_perfect_balanced_as_upper_bound = true;
    bool allow_overfitting = false;
    while (true) {
      HypernodeID hn = kInvalidHypernode;
      Gain gain = kInvalidGain;

      // The greedy initial partitioner has 3 different stages. In the first, we use the perfect
      // balanced part weight as upper bound for the block weights. Once we reach the block weight
      // limit, we release the upper bound and use the maximum allowed block weight as new upper bound.
      // Once we are not able to assign any vertex to a block, we allow overfitting, which effectively
      // allows to violate the balance constraint.
      if ( !PQSelectionPolicy::pop(hg, kway_pq, hn, to, gain, use_perfect_balanced_as_upper_bound) ) {
        if ( use_perfect_balanced_as_upper_bound ) {
          enableAllPQs(_context.partition.k, kway_pq);
          use_perfect_balanced_as_upper_bound = false;
          continue;
        } else if ( !allow_overfitting ) {
          enableAllPQs(_context.partition.k, kway_pq);
          allow_overfitting = true;
          continue;
        } else {
          break;
        }
      }

      ASSERT(hn != kInvalidHypernode);
      ASSERT(to != kInvalidPartition);
      ASSERT(hg.partID(hn) == kInvalidPartition);

      if ( allow_overfitting || fitsIntoBlock(hg, hn, to, use_perfect_balanced_as_upper_bound) ) {
        hg.setNodePart(hn, to);
        insertAndUpdateVerticesAfterMove(hg, kway_pq, hyperedges_in_queue, hn, to);
      } else {
        kway_pq.insert(hg.originalNodeID(hn), to, gain);
        kway_pq.disablePart(to);
      }
    }

    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(end - start).count();
    _ip_data.commit(_algorithm, time);
    return nullptr;
  }

 private:
  bool fitsIntoBlock(HyperGraph& hypergraph,
                     const HypernodeID hn,
                     const PartitionID block,
                     const bool use_perfect_balanced_as_upper_bound) const {
    ASSERT(block != kInvalidPartition && block < _context.partition.k);
    const HyperedgeWeight upper_bound = use_perfect_balanced_as_upper_bound ?
      _context.partition.perfect_balance_part_weights[block] : _context.partition.max_part_weights[block];
    return hypergraph.localPartWeight(block) + hypergraph.nodeWeight(hn) <=
      upper_bound;
  }

  void insertVertexIntoPQ(const HyperGraph& hypergraph,
                          KWayPriorityQueue& pq,
                          const HypernodeID hn,
                          const PartitionID to) {
    const HypernodeID original_id = hypergraph.originalNodeID(hn);
    ASSERT(to != kInvalidPartition && to < _context.partition.k);
    ASSERT(hypergraph.partID(hn) == kInvalidPartition);
    ASSERT(!pq.contains(original_id, to));

    const Gain gain = GainPolicy::calculateGain(hypergraph, hn, to);
    pq.insert(original_id, to, gain);
    if ( !pq.isEnabled(to) ) {
      pq.enablePart(to);
    }

    ASSERT(pq.contains(original_id, to));
    ASSERT(pq.isEnabled(to));
  }

  void insertUnassignedVertexIntoPQ(const HyperGraph& hypergraph,
                                    KWayPriorityQueue& pq,
                                    const PartitionID to) {
    const HypernodeID unassigned_hn = _ip_data.get_unassigned_hypernode();
    if ( unassigned_hn != kInvalidHypernode ) {
      insertVertexIntoPQ(hypergraph, pq, unassigned_hn, to);
    }
  }

  void insertAndUpdateVerticesAfterMove(const HyperGraph& hypergraph,
                                        KWayPriorityQueue& pq,
                                        kahypar::ds::FastResetFlagArray<>& hyperedges_in_queue,
                                        const HypernodeID hn,
                                        const PartitionID to) {
    ASSERT(to != kInvalidPartition && to < _context.partition.k);
    ASSERT(hypergraph.partID(hn) == to);

    // Perform delta gain updates
    GainPolicy::deltaGainUpdate(hypergraph, pq, hn, to);

    // Remove moved hypernode hn from all PQs
    const HypernodeID original_id = hypergraph.originalNodeID(hn);
    for ( PartitionID block = 0; block < hypergraph.k(); ++block ) {
      if ( pq.contains(original_id, block) ) {

        // Prevent that PQ becomes empty
        if ( to != block && pq.size(block) == 1 ) {
          insertUnassignedVertexIntoPQ(hypergraph, pq, block);
        }

        pq.remove(original_id, block);
      }
    }

    // Insert all adjacent hypernodes of the moved vertex into PQ of block to
    for ( const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
      const HyperedgeID original_he_id = hypergraph.originalEdgeID(he);
      if ( !hyperedges_in_queue[to * hypergraph.initialNumEdges() + original_he_id] ) {
        for ( const HypernodeID& pin : hypergraph.pins(he) ) {
          const HypernodeID original_pin_id = hypergraph.originalNodeID(pin);
          if ( hypergraph.partID(pin) == kInvalidPartition && !pq.contains(original_pin_id, to) ) {
            insertVertexIntoPQ(hypergraph, pq, pin, to);
          }
        }
        hyperedges_in_queue.set(to * hypergraph.initialNumEdges() + original_he_id, true);
      }
    }

    // Prevent that PQ becomes empty
    if ( pq.size(to) == 0 ) {
      insertUnassignedVertexIntoPQ(hypergraph, pq, to);
    }
  }

  void enableAllPQs(const PartitionID k, KWayPriorityQueue& pq) {
    for ( PartitionID block = 0; block < k; ++block ) {
      pq.enablePart(block);
    }
  }

  const InitialPartitioningAlgorithm _algorithm;
  InitialPartitioningDataContainer& _ip_data;
  const Context& _context;
};

template <typename TypeTraits, typename GainPolicy, typename PQSelectionPolicy>
PartitionID GreedyInitialPartitionerT<TypeTraits, GainPolicy, PQSelectionPolicy>::kInvalidPartition = -1;
template <typename TypeTraits, typename GainPolicy, typename PQSelectionPolicy>
HypernodeID GreedyInitialPartitionerT<TypeTraits, GainPolicy, PQSelectionPolicy>::kInvalidHypernode = std::numeric_limits<HypernodeID>::max();
template <typename TypeTraits, typename GainPolicy, typename PQSelectionPolicy>
Gain GreedyInitialPartitionerT<TypeTraits, GainPolicy, PQSelectionPolicy>::kInvalidGain = std::numeric_limits<Gain>::min();

template <typename GainPolicy, typename PQSelectionPolicy>
using GreedyInitialPartitioner = GreedyInitialPartitionerT<GlobalTypeTraits, GainPolicy, PQSelectionPolicy>;

} // namespace mt_kahypar

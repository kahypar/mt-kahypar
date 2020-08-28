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
#include "mt-kahypar/partition/initial_partitioning/flat/policies/pseudo_peripheral_start_nodes.h"

namespace mt_kahypar {
template<typename GainPolicy,
         typename PQSelectionPolicy>
class GreedyInitialPartitioner : public tbb::task {

  using DeltaFunction = std::function<void (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID)>;
  #define NOOP_FUNC [] (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID) { }

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

 public:
  GreedyInitialPartitioner(const InitialPartitioningAlgorithm algorithm,
                            InitialPartitioningDataContainer& ip_data,
                            const Context& context,
                            const int seed) :
    _algorithm(algorithm),
    _ip_data(ip_data),
    _context(context),
    _default_block(PQSelectionPolicy::getDefaultBlock()),
    _rng(seed) { }

  tbb::task* execute() override {
    if ( _ip_data.should_initial_partitioner_run(_algorithm) ) {
      HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
      PartitionedHypergraph& hg = _ip_data.local_partitioned_hypergraph();
      KWayPriorityQueue& kway_pq = _ip_data.local_kway_priority_queue();
      kahypar::ds::FastResetFlagArray<>& hyperedges_in_queue =
        _ip_data.local_hyperedge_fast_reset_flag_array();

      // Experiments have shown that some pq selection policies work better
      // if we preassign all vertices to a block and than execute the greedy
      // initial partitioner. E.g. the round-robin variant leaves the hypernode
      // unassigned, but the global and sequential strategy both preassign
      // all vertices to block 1 before initial partitioning.
      if ( _default_block != kInvalidPartition ) {
        ASSERT(_default_block < _context.partition.k);
        kway_pq.disablePart(_default_block);
        for ( const HypernodeID& hn : hg.nodes() ) {
          hg.setNodePart(hn, _default_block);
        }
      }

      // Insert start vertices into its corresponding PQs
      _ip_data.reset_unassigned_hypernodes();
      parallel::scalable_vector<HypernodeID> start_nodes =
        PseudoPeripheralStartNodes::computeStartNodes(_ip_data, _context, _default_block, _rng);
      ASSERT(static_cast<size_t>(_context.partition.k) == start_nodes.size());
      kway_pq.clear();
      for ( PartitionID block = 0; block < _context.partition.k; ++block ) {
        if ( block != _default_block ) {
          insertVertexIntoPQ(hg, kway_pq, start_nodes[block], block);
        }
      }

      hyperedges_in_queue.reset();
      PartitionID to = kInvalidPartition;
      bool use_perfect_balanced_as_upper_bound = true;
      bool allow_overfitting = false;
      while (true) {
        // If our default block has a weight less than the perfect balanced block weight
        // we terminate greedy initial partitioner in order to prevent that the default block
        // becomes underloaded.
        if ( _default_block != kInvalidPartition &&
            hg.partWeight(_default_block) <
            _context.partition.perfect_balance_part_weights[_default_block] ) {
          break;
        }

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
        ASSERT(to != _default_block);
        ASSERT(hg.partID(hn) == _default_block);

        if ( allow_overfitting || fitsIntoBlock(hg, hn, to, use_perfect_balanced_as_upper_bound) ) {
          if ( _default_block != kInvalidPartition ) {
            hg.changeNodePart(hn, _default_block, to, NOOP_FUNC);
          } else {
            hg.setNodePart(hn, to);
          }
          insertAndUpdateVerticesAfterMove(hg, kway_pq, hyperedges_in_queue, hn, _default_block, to);
        } else {
          kway_pq.insert(hn, to, gain);
          kway_pq.disablePart(to);
        }
      }

      HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
      double time = std::chrono::duration<double>(end - start).count();
      _ip_data.commit(_algorithm, time);
    }
    return nullptr;
  }

 private:
  bool fitsIntoBlock(PartitionedHypergraph& hypergraph,
                     const HypernodeID hn,
                     const PartitionID block,
                     const bool use_perfect_balanced_as_upper_bound) const {
    ASSERT(block != kInvalidPartition && block < _context.partition.k);
    const HyperedgeWeight upper_bound = use_perfect_balanced_as_upper_bound ?
      _context.partition.perfect_balance_part_weights[block] : _context.partition.max_part_weights[block];
    return hypergraph.partWeight(block) + hypergraph.nodeWeight(hn) <=
      upper_bound;
  }

  void insertVertexIntoPQ(const PartitionedHypergraph& hypergraph,
                          KWayPriorityQueue& pq,
                          const HypernodeID hn,
                          const PartitionID to) {
    ASSERT(to != kInvalidPartition && to < _context.partition.k);
    ASSERT(hypergraph.partID(hn) == _default_block, V(hypergraph.partID(hn)) << V(_default_block));
    ASSERT(!pq.contains(hn, to));

    const Gain gain = GainPolicy::calculateGain(hypergraph, hn, to);
    pq.insert(hn, to, gain);
    if ( !pq.isEnabled(to) ) {
      pq.enablePart(to);
    }

    ASSERT(pq.contains(hn, to));
    ASSERT(pq.isEnabled(to));
  }

  void insertUnassignedVertexIntoPQ(const PartitionedHypergraph& hypergraph,
                                    KWayPriorityQueue& pq,
                                    const PartitionID to) {
    ASSERT(to != _default_block);
    const HypernodeID unassigned_hn = _ip_data.get_unassigned_hypernode(_default_block);
    if ( unassigned_hn != kInvalidHypernode ) {
      insertVertexIntoPQ(hypergraph, pq, unassigned_hn, to);
    }
  }

  void insertAndUpdateVerticesAfterMove(const PartitionedHypergraph& hypergraph,
                                        KWayPriorityQueue& pq,
                                        kahypar::ds::FastResetFlagArray<>& hyperedges_in_queue,
                                        const HypernodeID hn,
                                        const PartitionID from,
                                        const PartitionID to) {
    ASSERT(to != kInvalidPartition && to < _context.partition.k);
    ASSERT(hypergraph.partID(hn) == to);

    // Perform delta gain updates
    GainPolicy::deltaGainUpdate(hypergraph, pq, hn, from, to);

    // Remove moved hypernode hn from all PQs
    for ( PartitionID block = 0; block < hypergraph.k(); ++block ) {
      if ( pq.contains(hn, block) ) {

        // Prevent that PQ becomes empty
        if ( to != block && pq.size(block) == 1 ) {
          insertUnassignedVertexIntoPQ(hypergraph, pq, block);
        }

        pq.remove(hn, block);
      }
    }

    // Insert all adjacent hypernodes of the moved vertex into PQ of block to
    for ( const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
      if ( !hyperedges_in_queue[to * hypergraph.initialNumEdges() + he] ) {
        for ( const HypernodeID& pin : hypergraph.pins(he) ) {
          if ( hypergraph.partID(pin) == _default_block && !pq.contains(pin, to) ) {
            insertVertexIntoPQ(hypergraph, pq, pin, to);
          }
        }
        hyperedges_in_queue.set(to * hypergraph.initialNumEdges() + he, true);
      }
    }

    // Prevent that PQ becomes empty
    if ( pq.size(to) == 0 ) {
      insertUnassignedVertexIntoPQ(hypergraph, pq, to);
    }
  }

  void enableAllPQs(const PartitionID k, KWayPriorityQueue& pq) {
    for ( PartitionID block = 0; block < k; ++block ) {
      if ( block != _default_block ) {
        pq.enablePart(block);
      }
    }
  }

  const InitialPartitioningAlgorithm _algorithm;
  InitialPartitioningDataContainer& _ip_data;
  const Context& _context;
  const PartitionID _default_block;
  std::mt19937 _rng;
};

} // namespace mt_kahypar

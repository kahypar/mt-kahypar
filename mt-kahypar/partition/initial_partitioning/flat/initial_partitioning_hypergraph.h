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

#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
template <typename TypeTraits>
class InitialPartitioningHypergraphT {
  using HyperGraph = typename TypeTraits::HyperGraph;
  using TBB = typename TypeTraits::TBB;

  static constexpr bool debug = false;
  static PartitionID kInvalidPartition;
  static HypernodeID kInvalidHypernode;

  struct PartitioningResult {
    explicit PartitioningResult(HyperedgeWeight objective,
                                double imbalance) :
      _objective(objective),
      _imbalance(imbalance) { }

    bool is_better(const PartitioningResult& other, const double epsilon) {
      bool improved_metric_within_balance = (other._imbalance <= epsilon &&
                                             other._objective < _objective);
      bool improved_balanced_of_imbalanced_partition = (_imbalance > epsilon &&
                                                        other._imbalance < _imbalance);
      bool equal_imbalance_but_better_metric = (_imbalance == other._imbalance &&
                                                other._objective < _objective);
      return improved_metric_within_balance ||
             improved_balanced_of_imbalanced_partition ||
             equal_imbalance_but_better_metric;
    }

    HyperedgeWeight _objective;
    double _imbalance;
  };

  struct LocalInitialPartitioningHypergraph {

    LocalInitialPartitioningHypergraph() :
      _hypergraph(),
      _mapping(),
      _partition(),
      _result(std::numeric_limits<HypernodeWeight>::max(), std::numeric_limits<double>::max()) { }

    LocalInitialPartitioningHypergraph(HyperGraph&& hypergraph,
                                       parallel::scalable_vector<HypernodeID>&& mapping) :
      _hypergraph(std::move(hypergraph)),
      _mapping(std::move(mapping)),
      _partition(_hypergraph.initialNumNodes(), kInvalidPartition),
      _result(std::numeric_limits<HypernodeWeight>::max(), std::numeric_limits<double>::max()) { }

    HyperGraph _hypergraph;
    parallel::scalable_vector<HypernodeID> _mapping;
    parallel::scalable_vector<PartitionID> _partition;
    PartitioningResult _result;
  };

  using ThreadLocalHypergraph = tbb::enumerable_thread_specific<LocalInitialPartitioningHypergraph>;
  using ThreadLocalUnassignedHypernodes = tbb::enumerable_thread_specific<parallel::scalable_vector<HypernodeID>>;

 public:
  InitialPartitioningHypergraphT(HyperGraph& hypergraph,
                                 const Context& context,
                                 const TaskGroupID& task_group_id) :
    _hg(hypergraph),
    _context(context),
    _task_group_id(task_group_id),
    _local_hg([&] {
      return copy_hypergraph();
    }),
    _local_hn_visited(_context.partition.k * _local_hg.local()._hypergraph.initialNumNodes()),
    _local_he_visited(_context.partition.k * _local_hg.local()._hypergraph.initialNumEdges()),
    _local_unassigned_hypernodes(),
    _local_unassigned_hypernode_pointer(std::numeric_limits<size_t>::max())  { }

  InitialPartitioningHypergraphT(const InitialPartitioningHypergraphT&) = delete;
  InitialPartitioningHypergraphT & operator= (const InitialPartitioningHypergraphT &) = delete;

  InitialPartitioningHypergraphT(InitialPartitioningHypergraphT&&) = delete;
  InitialPartitioningHypergraphT & operator= (InitialPartitioningHypergraphT &&) = delete;

  const HyperGraph& global_hypergraph() const {
    return _hg;
  }

  HyperGraph& local_hypergraph() {
    return _local_hg.local()._hypergraph;
  }

  kahypar::ds::FastResetFlagArray<>& local_hypernode_fast_reset_flag_array() {
    return _local_hn_visited.local();
  }

  kahypar::ds::FastResetFlagArray<>& local_hyperedge_fast_reset_flag_array() {
    return _local_he_visited.local();
  }

  void reset_unassigned_hypernodes() {
    parallel::scalable_vector<HypernodeID>& unassigned_hypernodes =
      _local_unassigned_hypernodes.local();
    size_t& unassigned_hypernode_pointer = _local_unassigned_hypernode_pointer.local();
    if ( unassigned_hypernode_pointer == std::numeric_limits<size_t>::max() ) {
      // In case the local unassigned hypernode vector was not initialized before
      // we initialize it here
      const HyperGraph& hypergraph = _local_hg.local()._hypergraph;
      for ( const HypernodeID& hn : hypergraph.nodes() ) {
        unassigned_hypernodes.push_back(hn);
      }
    }
    unassigned_hypernode_pointer = unassigned_hypernodes.size();
  }

  HypernodeID get_unassigned_hypernode() {
    const HyperGraph& hypergraph = _local_hg.local()._hypergraph;
    parallel::scalable_vector<HypernodeID>& unassigned_hypernodes =
      _local_unassigned_hypernodes.local();
    size_t& unassigned_hypernode_pointer = _local_unassigned_hypernode_pointer.local();
    ASSERT(unassigned_hypernodes.size() > 0);
    ASSERT(unassigned_hypernode_pointer <= unassigned_hypernodes.size());

    while ( unassigned_hypernode_pointer > 0 ) {
      const HypernodeID current_hn = unassigned_hypernodes[0];
      // In case the current hypernode is unassigned we return it
      if ( hypergraph.partID(current_hn) == kInvalidPartition ) {
        return current_hn;
      }
      // In case the hypernode on the first position is already assigned,
      // we swap it to end of the unassigned hypernode vector and decrement
      // the pointer such that we will not visit it again
      std::swap(unassigned_hypernodes[0], unassigned_hypernodes[--unassigned_hypernode_pointer]);
    }
    return kInvalidHypernode;
  }

  // ! Only for testing
  HypernodeID map_hypernode_to_local_hypergraph(const HypernodeID hn) {
    LocalInitialPartitioningHypergraph& local_hg = _local_hg.local();
    return local_hg._hypergraph.globalNodeID(local_hg._mapping[_hg.originalNodeID(hn)]);
  }

  /*!
   * Commits the current partition computed on the local hypergraph. Partition replaces
   * the best local partition, if it has a better quality (or better imbalance).
   * Partition on the local hypergraph is resetted afterwards.
   */
  void commit() {
    LocalInitialPartitioningHypergraph& ip_hypergraph = _local_hg.local();
    ip_hypergraph._hypergraph.initializeNumCutHyperedges();
    ip_hypergraph._hypergraph.updateGlobalPartInfos();

    PartitioningResult result(
      metrics::objective(ip_hypergraph._hypergraph, _context.partition.objective),
      metrics::imbalance(ip_hypergraph._hypergraph, _context));

    if ( ip_hypergraph._result.is_better(result, _context.partition.epsilon) ) {
      for ( const HypernodeID& hn : ip_hypergraph._hypergraph.nodes() ) {
        const HypernodeID original_id = ip_hypergraph._hypergraph.originalNodeID(hn);
        const PartitionID part_id = ip_hypergraph._hypergraph.partID(hn);
        ASSERT(original_id < ip_hypergraph._partition.size());
        ASSERT(part_id != kInvalidPartition);
        ip_hypergraph._partition[original_id] = part_id;
      }
      ip_hypergraph._result = std::move(result);
    }

    ip_hypergraph._hypergraph.resetPartition();
  }

  /*!
   * Determines the best partition computed by all threads and applies it to
   * the hypergraph. Note, this function is not thread-safe and should be called
   * if no other thread using that object operates on it.
   * Note, object is afterwards not usable any more.
   */
  void apply() {
    // Determine best partition
    LocalInitialPartitioningHypergraph best;
    for ( LocalInitialPartitioningHypergraph& partition : _local_hg ) {
      if ( best._result.is_better(partition._result, _context.partition.epsilon) ) {
        best = std::move(partition);
      }
    }

    // Applies best partition to hypergraph
    for ( const HypernodeID& hn : _hg.nodes() ) {
      const HypernodeID original_id = _hg.originalNodeID(hn);
      ASSERT(original_id < best._mapping.size());
      const PartitionID part_id = best._partition[best._mapping[original_id]];
      ASSERT(part_id != kInvalidPartition && part_id < _hg.k());
      _hg.setNodePart(hn, part_id);
    }
    _hg.initializeNumCutHyperedges();
    _hg.updateGlobalPartInfos();
    ASSERT(best._result._objective == metrics::objective(_hg, _context.partition.objective));
  }

 private:
  LocalInitialPartitioningHypergraph copy_hypergraph() {
    auto tmp_hg = _hg.copy(_context.partition.k, _task_group_id);
    return LocalInitialPartitioningHypergraph(std::move(tmp_hg.first), std::move(tmp_hg.second));
  }

  HyperGraph& _hg;
  const Context& _context;
  const TaskGroupID& _task_group_id;

  ThreadLocalHypergraph _local_hg;
  ThreadLocalFastResetFlagArray _local_hn_visited;
  ThreadLocalFastResetFlagArray _local_he_visited;
  ThreadLocalUnassignedHypernodes _local_unassigned_hypernodes;
  tbb::enumerable_thread_specific<size_t> _local_unassigned_hypernode_pointer;
};

template <typename TypeTraits>
PartitionID InitialPartitioningHypergraphT<TypeTraits>::kInvalidPartition = -1;
template <typename TypeTraits>
HypernodeID InitialPartitioningHypergraphT<TypeTraits>::kInvalidHypernode = std::numeric_limits<HypernodeID>::max();
} // namespace mt_kahypar
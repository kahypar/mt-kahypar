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

#include <sstream>

#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/label_propagation_refiner.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {

template <typename TypeTraits>
class InitialPartitioningDataContainerT {
  using HyperGraph = typename TypeTraits::HyperGraph;
  using TBB = typename TypeTraits::TBB;

  static constexpr bool debug = false;
  static PartitionID kInvalidPartition;
  static HypernodeID kInvalidHypernode;

  struct PartitioningResult {
    explicit PartitioningResult(InitialPartitioningAlgorithm algorithm,
                                HyperedgeWeight objective,
                                double imbalance) :
      _algorithm(algorithm),
      _objective(objective),
      _imbalance(imbalance) { }

    bool is_other_better(const PartitioningResult& other, const double epsilon) {
      bool equal_metric = other._objective == _objective;
      bool improved_metric = other._objective < _objective;
      bool improved_imbalance = other._imbalance < _imbalance;
      bool is_feasible = _imbalance <= epsilon;
      bool is_other_feasible = other._imbalance <= epsilon;
      return ( improved_metric && (is_other_feasible || improved_imbalance) ) ||
             ( equal_metric && improved_imbalance ) ||
             ( is_other_feasible && !is_feasible ) ||
             ( improved_metric && !is_other_feasible && !is_feasible );
    }

    std::string str() const {
      std::stringstream ss;
      ss << "Algorithm = " << _algorithm << ", "
         << "Objective = " << _objective << ", "
         << "Imbalance = " << _imbalance;
      return ss.str();
    }

    InitialPartitioningAlgorithm _algorithm;
    HyperedgeWeight _objective;
    double _imbalance;
  };

  struct LocalInitialPartitioningHypergraph {
    using LabelPropagationKm1Refiner = LabelPropagationRefinerT<TypeTraits, NoExecutionPolicy, Km1Policy>;
    using LabelPropagationCutRefiner = LabelPropagationRefinerT<TypeTraits, NoExecutionPolicy, CutPolicy>;

    LocalInitialPartitioningHypergraph(HyperGraph&& hypergraph,
                                       const Context& context,
                                       parallel::scalable_vector<HypernodeID>&& mapping,
                                       const TaskGroupID task_group_id) :
      _hypergraph(std::move(hypergraph)),
      _context(context),
      _mapping(std::move(mapping)),
      _partition(_hypergraph.initialNumNodes(), kInvalidPartition),
      _result(InitialPartitioningAlgorithm::UNDEFINED,
              std::numeric_limits<HypernodeWeight>::max(),
              std::numeric_limits<double>::max()),
      _label_propagation(nullptr) {

      if ( _context.refinement.label_propagation.algorithm != LabelPropagationAlgorithm::do_nothing ) {
        if ( _context.partition.objective == kahypar::Objective::km1 ) {
          _label_propagation = std::make_unique<LabelPropagationKm1Refiner>(
            _hypergraph, _context, task_group_id);
        } else if ( _context.partition.objective == kahypar::Objective::cut ) {
          _label_propagation = std::make_unique<LabelPropagationCutRefiner>(
            _hypergraph, _context, task_group_id);
        }
      }
    }

    void commit(const InitialPartitioningAlgorithm algorithm) {
      ASSERT([&]() {
          for (const HypernodeID& hn : _hypergraph.nodes()) {
            if (_hypergraph.partID(hn) == kInvalidPartition) {
              return false;
            }
          }
          return true;
        } (), "There are unassigned hypernodes!");

      _hypergraph.initializeNumCutHyperedges();
      _hypergraph.updateGlobalPartInfos();

      kahypar::Metrics current_metric = {
        metrics::hyperedgeCut(_hypergraph), metrics::km1(_hypergraph),
        metrics::imbalance(_hypergraph, _context) };

      if ( _label_propagation ) {
        _label_propagation->refine({}, current_metric);
      }

      PartitioningResult result(algorithm,
        current_metric.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
        current_metric.imbalance);

      DBG << "Thread ID:" << std::this_thread::get_id()
          << "- CPU ID:" << sched_getcpu()
          << "[" << result.str() << "]";

      if ( _result.is_other_better(result, _context.partition.epsilon) ) {
        for ( const HypernodeID& hn : _hypergraph.nodes() ) {
          const HypernodeID original_id = _hypergraph.originalNodeID(hn);
          const PartitionID part_id = _hypergraph.partID(hn);
          ASSERT(original_id < _partition.size());
          ASSERT(part_id != kInvalidPartition);
          _partition[original_id] = part_id;
        }
        _result = std::move(result);
      }

      _hypergraph.resetPartition();
    }

    HyperGraph _hypergraph;
    const Context& _context;
    parallel::scalable_vector<HypernodeID> _mapping;
    parallel::scalable_vector<PartitionID> _partition;
    PartitioningResult _result;
    std::unique_ptr<IRefiner> _label_propagation;
  };

  using ThreadLocalHypergraph = tbb::enumerable_thread_specific<LocalInitialPartitioningHypergraph>;
  using ThreadLocalUnassignedHypernodes = tbb::enumerable_thread_specific<parallel::scalable_vector<HypernodeID>>;

 public:
  InitialPartitioningDataContainerT(HyperGraph& hypergraph,
                                    const Context& context,
                                    const TaskGroupID task_group_id) :
    _hg(hypergraph),
    _context(context),
    _task_group_id(task_group_id),
    _local_hg([&] {
      return copy_hypergraph();
    }),
    _local_hn_visited(_context.partition.k * _local_hg.local()._hypergraph.initialNumNodes()),
    _local_he_visited(_context.partition.k * _local_hg.local()._hypergraph.initialNumEdges()),
    _local_unassigned_hypernodes(),
    _local_unassigned_hypernode_pointer(std::numeric_limits<size_t>::max())  {
    // Setup Label Propagation Refiner Config for Initial Partitioning
    _context.refinement.label_propagation.part_weight_update_frequency = std::numeric_limits<size_t>::max();
    _context.refinement.label_propagation.numa_aware = false;
    _context.refinement.label_propagation.localized = false;
    _context.refinement.label_propagation.execution_policy = ExecutionType::none;
    _context.refinement.label_propagation.execute_always = true;
    _context.refinement.label_propagation.execute_sequential = true;
  }

  InitialPartitioningDataContainerT(const InitialPartitioningDataContainerT&) = delete;
  InitialPartitioningDataContainerT & operator= (const InitialPartitioningDataContainerT &) = delete;

  InitialPartitioningDataContainerT(InitialPartitioningDataContainerT&&) = delete;
  InitialPartitioningDataContainerT & operator= (InitialPartitioningDataContainerT &&) = delete;

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
  void commit(const InitialPartitioningAlgorithm algorithm) {
    _local_hg.local().commit(algorithm);
  }

  /*!
   * Determines the best partition computed by all threads and applies it to
   * the hypergraph. Note, this function is not thread-safe and should be called
   * if no other thread using that object operates on it.
   */
  void apply() {
    // Determine best partition
    LocalInitialPartitioningHypergraph* best = nullptr;
    LocalInitialPartitioningHypergraph* worst = nullptr;
    LocalInitialPartitioningHypergraph* best_imbalance = nullptr;
    LocalInitialPartitioningHypergraph* best_objective = nullptr;
    for ( LocalInitialPartitioningHypergraph& partition : _local_hg ) {
      if ( !best || best->_result.is_other_better(partition._result, _context.partition.epsilon) ) {
        best = &partition;
      }
      if ( !worst || !worst->_result.is_other_better(partition._result, _context.partition.epsilon) ) {
        worst = &partition;
      }
      if ( !best_imbalance || best_imbalance->_result._imbalance > partition._result._imbalance ||
           (best_imbalance->_result._imbalance == partition._result._imbalance &&
            best_objective->_result._objective > partition._result._objective)) {
        best_imbalance = &partition;
      }
      if ( !best_objective || best_objective->_result._objective > partition._result._objective ) {
        best_objective = &partition;
      }
    }

    ASSERT(best);
    ASSERT(worst);
    ASSERT(best_imbalance);
    ASSERT(best_objective);
    DBG << "Num Vertices =" << _hg.initialNumNodes() << ", Num Edges =" << _hg.initialNumEdges()
        << ", k =" << _context.partition.k << ", epsilon =" << _context.partition.epsilon;
    DBG << "Best Partition                [" << best->_result.str() << "]";
    DBG << "Worst Partition               [" << worst->_result.str() << "]";
    DBG << "Best Balanced Partition       [" << best_imbalance->_result.str() << "]";
    DBG << "Partition with Best Objective [" << best_objective->_result.str() << "]";

    // Applies best partition to hypergraph
    for ( const HypernodeID& hn : _hg.nodes() ) {
      const HypernodeID original_id = _hg.originalNodeID(hn);
      ASSERT(original_id < best->_mapping.size());
      const PartitionID part_id = best->_partition[best->_mapping[original_id]];
      ASSERT(part_id != kInvalidPartition && part_id < _hg.k());
      _hg.setNodePart(hn, part_id);
    }

    _hg.initializeNumCutHyperedges();
    _hg.updateGlobalPartInfos();
    ASSERT(best->_result._objective == metrics::objective(_hg, _context.partition.objective),
           V(best->_result._objective) << V(metrics::objective(_hg, _context.partition.objective)));
  }

 private:
  LocalInitialPartitioningHypergraph copy_hypergraph() {
    auto tmp_hg = _hg.copy_sequential(_context.partition.k);
    return LocalInitialPartitioningHypergraph(
      std::move(tmp_hg.first), _context, std::move(tmp_hg.second), _task_group_id);
  }

  HyperGraph& _hg;
  Context _context;
  const TaskGroupID _task_group_id;

  ThreadLocalHypergraph _local_hg;
  ThreadLocalFastResetFlagArray _local_hn_visited;
  ThreadLocalFastResetFlagArray _local_he_visited;
  ThreadLocalUnassignedHypernodes _local_unassigned_hypernodes;
  tbb::enumerable_thread_specific<size_t> _local_unassigned_hypernode_pointer;
};

template <typename TypeTraits>
PartitionID InitialPartitioningDataContainerT<TypeTraits>::kInvalidPartition = -1;
template <typename TypeTraits>
HypernodeID InitialPartitioningDataContainerT<TypeTraits>::kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

using InitialPartitioningDataContainer = InitialPartitioningDataContainerT<GlobalTypeTraits>;
} // namespace mt_kahypar
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
#include <mutex>

#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/partition/initial_partitioning/flat/initial_partitioning_commons.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/initial_partitioning_stats.h"
#include "mt-kahypar/partition/refinement/fm/sequential_twoway_fm_refiner.h"


namespace mt_kahypar {

class InitialPartitioningDataContainer {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  // ! Contains information about the best thread local partition
  struct PartitioningResult {
    explicit PartitioningResult(InitialPartitioningAlgorithm algorithm,
                                HyperedgeWeight objective_ip,
                                HyperedgeWeight objective,
                                double imbalance) :
      _algorithm(algorithm),
      _objective_ip(objective_ip),
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
             ( improved_imbalance && !is_other_feasible && !is_feasible );
    }

    std::string str() const {
      std::stringstream ss;
      ss << "Algorithm = " << _algorithm << ", "
         << "Objective IP = " << _objective_ip << ", "
         << "Objective = " << _objective << ", "
         << "Imbalance = " << _imbalance;
      return ss.str();
    }

    InitialPartitioningAlgorithm _algorithm;
    HyperedgeWeight _objective_ip;
    HyperedgeWeight _objective;
    double _imbalance;
  };

  // ! Aggregates global stats about the partitions produced by an specific
  // ! initial partitioning algorithm.
  struct InitialPartitioningRunStats {
    explicit InitialPartitioningRunStats(InitialPartitioningAlgorithm algo) :
      algorithm(algo),
      average_quality(0.0),
      sum_of_squares(0.0),
      n(0),
      best_quality(std::numeric_limits<HyperedgeWeight>::max()) { }

    void add_run(const HyperedgeWeight quality) {
      ++n;
      // Incremental update standard deviation
      // Incremental update average quality
      const double old_average_quality = average_quality;
      average_quality += static_cast<double>(quality - average_quality) / n;
      sum_of_squares += (quality - old_average_quality) * (quality - average_quality);
      if ( quality < best_quality ) {
        best_quality = quality;
      }
    }

    double stddev() const {
      return n == 1 ? 0 : std::sqrt(sum_of_squares / ( n - 1 ));
    }

    InitialPartitioningAlgorithm algorithm;
    double average_quality;
    long double sum_of_squares;
    size_t n;
    HyperedgeWeight best_quality;
  };

  // ! Aggregates global stats of all initial partitioning algorithms.
  // ! Additionally it provides a function that decides whether it is
  // ! beneficial to perform additional runs of a specific initial
  // ! partitioning algorithm based on its previous runs (see
  // ! should_initial_partitioner_run(...)).
  struct GlobalInitialPartitioningStats {

    explicit GlobalInitialPartitioningStats(const Context& context) :
      _stat_mutex(),
      _context(context),
      _stats(),
      _best_quality(std::numeric_limits<HyperedgeWeight>::max()) {
      const uint8_t num_initial_partitioner = static_cast<uint8_t>(InitialPartitioningAlgorithm::UNDEFINED);
      for ( uint8_t algo = 0; algo < num_initial_partitioner; ++algo ) {
        _stats.emplace_back(static_cast<InitialPartitioningAlgorithm>(algo));
      }
    }

    void add_run(const InitialPartitioningAlgorithm algorithm,
                 const HyperedgeWeight quality,
                 const bool is_feasible) {
      std::lock_guard<std::mutex> _lock(_stat_mutex);
      const uint8_t algo_idx = static_cast<uint8_t>(algorithm);
      _stats[algo_idx].add_run(quality);
      if ( is_feasible && quality < _best_quality ) {
        _best_quality = quality;
      }
    }

    // ! Decides whether it is beneficial to perform further runs of a specific
    // ! initial partitioning algorithm. Function assumes that the quality produced
    // ! by a partitioner follows a normal distribution. In that case, approx. 95%
    // ! of the partitions produced by an initial partitioner have a quality between
    // ! avg_quality - 2 * stddev_quality and avg_quality + 2 * stddev_quality. If
    // ! avg_quality - 2 * stddev_quality is greater than the best partition produced
    // ! so far then we say that the probability that the corresponding initial
    // ! partitioner produce a new global best partition is too low and prohibit further
    // ! runs of that partitioner.
    bool should_initial_partitioner_run(const InitialPartitioningAlgorithm algorithm) const {
      const uint8_t algo_idx = static_cast<uint8_t>(algorithm);
      return !_context.initial_partitioning.use_adaptive_ip_runs ||
             _stats[algo_idx].n < _context.initial_partitioning.min_adaptive_ip_runs ||
             _stats[algo_idx].average_quality - 2.0 * _stats[algo_idx].stddev() <= _best_quality;
    }

    std::mutex _stat_mutex;
    const Context& _context;
    parallel::scalable_vector<InitialPartitioningRunStats> _stats;
    HyperedgeWeight _best_quality;
  };

  struct LocalInitialPartitioningHypergraph {

    LocalInitialPartitioningHypergraph(Hypergraph& hypergraph,
                                       const Context& context,
                                       const TaskGroupID task_group_id,
                                       GlobalInitialPartitioningStats& global_stats,
                                       const bool disable_fm) :
      _partitioned_hypergraph(context.partition.k, hypergraph),
      _context(context),
      _global_stats(global_stats),
      _partition(hypergraph.initialNumNodes(), kInvalidPartition),
      _result(InitialPartitioningAlgorithm::UNDEFINED,
              std::numeric_limits<HypernodeWeight>::max(),
              std::numeric_limits<HypernodeWeight>::max(),
              std::numeric_limits<double>::max()),
      _label_propagation(nullptr),
      _twoway_fm(nullptr),
      _stats() {

      for ( uint8_t algo = 0; algo < static_cast<size_t>(InitialPartitioningAlgorithm::UNDEFINED); ++algo ) {
        _stats.emplace_back(static_cast<InitialPartitioningAlgorithm>(algo));
      }

      if ( _context.partition.k == 2 && !disable_fm ) {
        // In case of a bisection we instantiate the 2-way FM refiner
        _twoway_fm = std::make_unique<SequentialTwoWayFmRefiner>(_partitioned_hypergraph, _context);
      } else if ( _context.refinement.label_propagation.algorithm != LabelPropagationAlgorithm::do_nothing ) {
        // In case of a direct-kway initial partition we instantiate the LP refiner
        _label_propagation = LabelPropagationFactory::getInstance().createObject(
          _context.refinement.label_propagation.algorithm, hypergraph, _context, task_group_id);
      }
    }

    void commit(const InitialPartitioningAlgorithm algorithm,
                const double time = 0.0) {
      ASSERT([&]() {
          for (const HypernodeID& hn : _partitioned_hypergraph.nodes()) {
            if (_partitioned_hypergraph.partID(hn) == kInvalidPartition) {
              return false;
            }
          }
          return true;
        } (), "There are unassigned hypernodes!");


      kahypar::Metrics current_metric = {
        metrics::hyperedgeCut(_partitioned_hypergraph, false),
        metrics::km1(_partitioned_hypergraph, false),
        metrics::imbalance(_partitioned_hypergraph, _context) };
      const HyperedgeWeight quality_before_refinement =
        current_metric.getMetric(kahypar::Mode::direct_kway, _context.partition.objective);

      refineCurrentPartition(current_metric);

      PartitioningResult result(algorithm, quality_before_refinement,
        current_metric.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
        current_metric.imbalance);

      DBG << "[" << result.str() << "]";

      // Aggregate Stats
      uint8_t algorithm_index = static_cast<uint8_t>(algorithm);
      _stats[algorithm_index].total_sum_quality += result._objective;
      _stats[algorithm_index].total_time += time;
      ++_stats[algorithm_index].total_calls;

      if ( _result.is_other_better(result, _context.partition.epsilon) ) {
        for ( const HypernodeID& hn : _partitioned_hypergraph.nodes() ) {
          const PartitionID part_id = _partitioned_hypergraph.partID(hn);
          ASSERT(hn < _partition.size());
          ASSERT(part_id != kInvalidPartition);
          _partition[hn] = part_id;
        }
        _result = std::move(result);
      }

      _global_stats.add_run(algorithm, current_metric.getMetric(kahypar::Mode::direct_kway,
        _context.partition.objective), current_metric.imbalance <= _context.partition.epsilon);
      _partitioned_hypergraph.resetPartition();
    }

    void performRefinementOnBestPartition() {
      kahypar::Metrics current_metric = {
        _result._objective,
        _result._objective,
        _result._imbalance };

      // Apply best partition to hypergraph
      for ( const HypernodeID& hn : _partitioned_hypergraph.nodes() ) {
        ASSERT(hn < _partition.size());
        ASSERT(_partitioned_hypergraph.partID(hn) == kInvalidPartition);
        _partitioned_hypergraph.setNodePart(hn, _partition[hn]);
      }

      HEAVY_INITIAL_PARTITIONING_ASSERT(
        current_metric.getMetric(kahypar::Mode::direct_kway, _context.partition.objective) ==
        metrics::objective(_partitioned_hypergraph, _context.partition.objective, false));

      refineCurrentPartition(current_metric);

      // Compare current best partition with refined partition
      PartitioningResult result(_result._algorithm,
        current_metric.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
        current_metric.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
        current_metric.imbalance);

      if ( _result.is_other_better(result, _context.partition.epsilon) ) {
        for ( const HypernodeID& hn : _partitioned_hypergraph.nodes() ) {
          const PartitionID part_id = _partitioned_hypergraph.partID(hn);
          ASSERT(hn < _partition.size());
          ASSERT(part_id != kInvalidPartition);
          _partition[hn] = part_id;
        }
        _result = std::move(result);
      }
    }

    void refineCurrentPartition(kahypar::Metrics& current_metric) {
      if ( _context.partition.k == 2 && _twoway_fm ) {
        bool improvement = true;
        for ( size_t i = 0; i < _context.initial_partitioning.fm_refinment_rounds && improvement; ++i ) {
          improvement = _twoway_fm->refine(current_metric);
        }
      } else if ( _label_propagation ) {
        _label_propagation->initialize(_partitioned_hypergraph);
        _label_propagation->refine(_partitioned_hypergraph, {},
          current_metric, std::numeric_limits<double>::max());
      }

      HEAVY_INITIAL_PARTITIONING_ASSERT(
        current_metric.getMetric(kahypar::Mode::direct_kway, _context.partition.objective) ==
        metrics::objective(_partitioned_hypergraph, _context.partition.objective, false));
    }

    void aggregate_stats(parallel::scalable_vector<utils::InitialPartitionerSummary>& main_stats) const {
      ASSERT(main_stats.size() == _stats.size());
      for ( size_t i = 0; i < _stats.size(); ++i ) {
        main_stats[i].add(_stats[i]);
      }
    }

    void freeInternalData() {
      tbb::parallel_invoke([&] {
        _partitioned_hypergraph.freeInternalData();
      }, [&] {
        parallel::free(_partition);
      });
    }

    PartitionedHypergraph _partitioned_hypergraph;
    const Context& _context;
    GlobalInitialPartitioningStats& _global_stats;
    parallel::scalable_vector<PartitionID> _partition;
    PartitioningResult _result;
    std::unique_ptr<IRefiner> _label_propagation;
    std::unique_ptr<SequentialTwoWayFmRefiner> _twoway_fm;
    parallel::scalable_vector<utils::InitialPartitionerSummary> _stats;
  };

  using ThreadLocalHypergraph = tbb::enumerable_thread_specific<LocalInitialPartitioningHypergraph>;
  using ThreadLocalUnassignedHypernodes = tbb::enumerable_thread_specific<parallel::scalable_vector<HypernodeID>>;

 public:
  InitialPartitioningDataContainer(PartitionedHypergraph& hypergraph,
                                    const Context& context,
                                    const TaskGroupID task_group_id,
                                    const bool disable_fm = false) :
    _partitioned_hg(hypergraph),
    _context(context),
    _task_group_id(task_group_id),
    _disable_fm(disable_fm),
    _global_stats(context),
    _local_hg([&] {
      return construct_local_partitioned_hypergraph();
    }),
    _local_kway_pq(_context.partition.k),
    _is_local_pq_initialized(false),
    _local_hn_visited(_context.partition.k * hypergraph.initialNumNodes()),
    _local_he_visited(_context.partition.k * hypergraph.initialNumEdges()),
    _local_unassigned_hypernodes(),
    _local_unassigned_hypernode_pointer(std::numeric_limits<size_t>::max())  {
    // Setup Label Propagation IRefiner Config for Initial Partitioning
    _context.refinement = _context.initial_partitioning.refinement;
    _context.refinement.label_propagation.execute_sequential = true;
  }

  InitialPartitioningDataContainer(const InitialPartitioningDataContainer&) = delete;
  InitialPartitioningDataContainer & operator= (const InitialPartitioningDataContainer &) = delete;

  InitialPartitioningDataContainer(InitialPartitioningDataContainer&&) = delete;
  InitialPartitioningDataContainer & operator= (InitialPartitioningDataContainer &&) = delete;

  ~InitialPartitioningDataContainer() {
    tbb::parallel_invoke([&] {
      parallel::parallel_free_thread_local_internal_data(
        _local_hg, [&](LocalInitialPartitioningHypergraph& local_hg) {
          local_hg.freeInternalData();
        });
    }, [&] {
      parallel::parallel_free_thread_local_internal_data(
        _local_unassigned_hypernodes, [&](parallel::scalable_vector<HypernodeID>& array) {
          parallel::free(array);
        });
    });
  }

  PartitionedHypergraph& global_partitioned_hypergraph() {
    return _partitioned_hg;
  }

  PartitionedHypergraph& local_partitioned_hypergraph() {
    return _local_hg.local()._partitioned_hypergraph;
  }

  KWayPriorityQueue& local_kway_priority_queue() {
    bool& is_local_pq_initialized = _is_local_pq_initialized.local();
    KWayPriorityQueue& local_kway_pq = _local_kway_pq.local();
    if ( !is_local_pq_initialized ) {
      local_kway_pq.initialize(local_partitioned_hypergraph().initialNumNodes());
      is_local_pq_initialized = true;
    }
    return local_kway_pq;
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
      const PartitionedHypergraph& hypergraph = local_partitioned_hypergraph();
      for ( const HypernodeID& hn : hypergraph.nodes() ) {
        unassigned_hypernodes.push_back(hn);
      }
      utils::Randomize::instance().shuffleVector(unassigned_hypernodes, sched_getcpu());
    }
    unassigned_hypernode_pointer = unassigned_hypernodes.size();
  }

  HypernodeID get_unassigned_hypernode(const PartitionID unassigned_block = kInvalidPartition) {
    const PartitionedHypergraph& hypergraph = local_partitioned_hypergraph();
    parallel::scalable_vector<HypernodeID>& unassigned_hypernodes =
      _local_unassigned_hypernodes.local();
    size_t& unassigned_hypernode_pointer = _local_unassigned_hypernode_pointer.local();
    ASSERT(unassigned_hypernodes.size() > 0);
    ASSERT(unassigned_hypernode_pointer <= unassigned_hypernodes.size());

    while ( unassigned_hypernode_pointer > 0 ) {
      const HypernodeID current_hn = unassigned_hypernodes[0];
      // In case the current hypernode is unassigned we return it
      if ( hypergraph.partID(current_hn) == unassigned_block ) {
        return current_hn;
      }
      // In case the hypernode on the first position is already assigned,
      // we swap it to end of the unassigned hypernode vector and decrement
      // the pointer such that we will not visit it again
      std::swap(unassigned_hypernodes[0], unassigned_hypernodes[--unassigned_hypernode_pointer]);
    }
    return kInvalidHypernode;
  }

  bool should_initial_partitioner_run(const InitialPartitioningAlgorithm algorithm) {
    return _global_stats.should_initial_partitioner_run(algorithm);
  }

  /*!
   * Commits the current partition computed on the local hypergraph. Partition replaces
   * the best local partition, if it has a better quality (or better imbalance).
   * Partition on the local hypergraph is resetted afterwards.
   */
  void commit(const InitialPartitioningAlgorithm algorithm, const double time = 0.0) {
    _local_hg.local().commit(algorithm, time);
  }

  /*!
   * Determines the best partition computed by all threads and applies it to
   * the hypergraph. Note, this function is not thread-safe and should be called
   * if no other thread using that object operates on it.
   */
  void apply() {
    // Initialize Stats
    parallel::scalable_vector<utils::InitialPartitionerSummary> stats;
    size_t number_of_threads = 0;
    for ( uint8_t algo = 0; algo < static_cast<size_t>(InitialPartitioningAlgorithm::UNDEFINED); ++algo ) {
      stats.emplace_back(static_cast<InitialPartitioningAlgorithm>(algo));
    }

    // Perform FM refinement on the best partition of each thread
    if ( _context.initial_partitioning.perform_refinement_on_best_partitions ) {
      tbb::task_group fm_refinement_group;
      for ( LocalInitialPartitioningHypergraph& partition : _local_hg ) {
        fm_refinement_group.run([&] {
          partition.performRefinementOnBestPartition();
        });
      }
      fm_refinement_group.wait();
    }

    // Determine best partition
    LocalInitialPartitioningHypergraph* best = nullptr;
    LocalInitialPartitioningHypergraph* worst = nullptr;
    LocalInitialPartitioningHypergraph* best_imbalance = nullptr;
    LocalInitialPartitioningHypergraph* best_objective = nullptr;
    for ( LocalInitialPartitioningHypergraph& partition : _local_hg ) {
      ++number_of_threads;
      partition.aggregate_stats(stats);
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
    DBG << "Num Vertices =" << _partitioned_hg.initialNumNodes()
        << ", Num Edges =" << _partitioned_hg.initialNumEdges()
        << ", k =" << _context.partition.k << ", epsilon =" << _context.partition.epsilon;
    DBG << "Best Partition                [" << best->_result.str() << "]";
    DBG << "Worst Partition               [" << worst->_result.str() << "]";
    DBG << "Best Balanced Partition       [" << best_imbalance->_result.str() << "]";
    DBG << "Partition with Best Objective [" << best_objective->_result.str() << "]";

    // Applies best partition to hypergraph
    _partitioned_hg.doParallelForAllNodes([&](const HypernodeID hn) {
      ASSERT(hn < best->_partition.size());
      const PartitionID part_id = best->_partition[hn];
      ASSERT(part_id != kInvalidPartition && part_id < _partitioned_hg.k());
      _partitioned_hg.setOnlyNodePart(hn, part_id);
    });
    _partitioned_hg.initializePartition(_task_group_id);

    utils::InitialPartitioningStats::instance().add_initial_partitioning_result(
      best->_result._algorithm, number_of_threads, stats);
    ASSERT(best->_result._objective == metrics::objective(_partitioned_hg, _context.partition.objective, false),
           V(best->_result._objective) << V(metrics::objective(_partitioned_hg, _context.partition.objective, false)));
  }

 private:
  LocalInitialPartitioningHypergraph construct_local_partitioned_hypergraph() {
    return LocalInitialPartitioningHypergraph(
      _partitioned_hg.hypergraph(), _context, _task_group_id, _global_stats, _disable_fm);
  }

  PartitionedHypergraph& _partitioned_hg;
  Context _context;
  const TaskGroupID _task_group_id;
  const bool _disable_fm;

  GlobalInitialPartitioningStats _global_stats;
  ThreadLocalHypergraph _local_hg;
  ThreadLocalKWayPriorityQueue _local_kway_pq;
  tbb::enumerable_thread_specific<bool> _is_local_pq_initialized;
  ThreadLocalFastResetFlagArray _local_hn_visited;
  ThreadLocalFastResetFlagArray _local_he_visited;
  ThreadLocalUnassignedHypernodes _local_unassigned_hypernodes;
  tbb::enumerable_thread_specific<size_t> _local_unassigned_hypernode_pointer;
};

} // namespace mt_kahypar
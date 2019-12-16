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

#include <libkahypar.h>

#include "kahypar/partition/context.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/kahypar.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {
template <typename TypeTraits>
class DirectInitialPartitionerT : public IInitialPartitioner {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;

  static constexpr bool debug = false;
  static constexpr bool kahypar_debug = false;
  static PartitionID kInvalidPartition;
  static HypernodeID kInvalidHypernode;

 public:
  DirectInitialPartitionerT(HyperGraph& hypergraph, const Context& context, const bool) :
    _hg(hypergraph),
    _context(context) { }

  DirectInitialPartitionerT(const DirectInitialPartitionerT&) = delete;
  DirectInitialPartitionerT(DirectInitialPartitionerT&&) = delete;
  DirectInitialPartitionerT & operator= (const DirectInitialPartitionerT &) = delete;
  DirectInitialPartitionerT & operator= (DirectInitialPartitionerT &&) = delete;

 private:
  void initialPartitionImpl() override final {
    kahypar_context_t* context = setupContext(_context, kahypar_debug);

    // Setup number of runs per thread
    std::vector<size_t> ip_runs;
    if (_context.initial_partitioning.runs <= _context.shared_memory.num_threads ||
        _context.initial_partitioning.call_kahypar_multiple_times) {
      // In case of call_kahypar_multiple_times is true, we make runs calls to KaHyPar
      // and call the initial partitioner of KaHyPar only ONE time ...
      ip_runs.assign(_context.initial_partitioning.runs, 1);
    } else {
      // otherwise, we call KaHyPar num_threads times and split the initial partitioning runs of
      // KaHyPar evenly among the threads.
      size_t runs_per_thread = _context.initial_partitioning.runs / _context.shared_memory.num_threads;
      ip_runs.assign(_context.shared_memory.num_threads, runs_per_thread);
      for (size_t i = 0; i < _context.initial_partitioning.runs % _context.shared_memory.num_threads; ++i) {
        ++ip_runs[i];
      }
    }
    ASSERT([&] {
        size_t runs = 0;
        for (const size_t& n : ip_runs) {
          runs += n;
        }
        return runs == _context.initial_partitioning.runs;
      } (), "Number of runs per thread does not sum up to total number of ip runs");

    // Setup node mapping to a continous range
    HypernodeID num_vertices = 0;
    std::vector<HypernodeID> node_mapping(_hg.initialNumNodes(), kInvalidHypernode);
    std::vector<HypernodeID> reverse_mapping;
    for (const HypernodeID& hn : _hg.nodes()) {
      ASSERT(_hg.originalNodeID(hn) < _hg.initialNumNodes());
      node_mapping[_hg.originalNodeID(hn)] = num_vertices++;
      reverse_mapping.emplace_back(hn);
    }

    // Split calls to initial partitioner evenly across numa nodes
    size_t num_ip_calls = ip_runs.size();
    std::vector<KaHyParPartitioningResult> results(num_ip_calls);
    KaHyParHypergraph kahypar_hypergraph = convertToKaHyParHypergraph(_hg, node_mapping);

    size_t used_numa_nodes = TBB::instance().num_used_numa_nodes();
    std::vector<int> numa_cpus_prefix_sum(used_numa_nodes + 1, 0);
    for (size_t i = 1; i <= used_numa_nodes; ++i) {
      numa_cpus_prefix_sum[i] = numa_cpus_prefix_sum[i - 1] +
                                TBB::instance().number_of_threads_on_numa_node(i - 1);
    }
    // Returns a numa node with probability of the number of active threads on that
    // numa node divided by the total number of active threads
    // => this ensures that the initial partitioning runs are evenly scattered across
    //    the numa nodes
    auto get_numa_node_for_execution = [&]() {
      int numa_node = utils::Randomize::instance().getRandomInt(0,
      numa_cpus_prefix_sum[used_numa_nodes - 1], sched_getcpu());
      for (size_t j = 0; j < used_numa_nodes; ++j) {
        if (numa_cpus_prefix_sum[j] <= numa_node && numa_node < numa_cpus_prefix_sum[j + 1]) {
          numa_node = j;
        }
      }
      ASSERT(numa_node < (int)used_numa_nodes);
      return numa_node;
    };

    for (size_t i = 0; i < num_ip_calls; ++i) {
      int numa_node = get_numa_node_for_execution();
      TBB::instance().numa_task_arena(numa_node).execute([&, i] {
          TBB::instance().numa_task_group(numa_node).run([&, i] {
            size_t seed = _context.partition.seed + i * _context.initial_partitioning.runs;
            results[i] = partitionWithKaHyPar(kahypar_hypergraph, context, ip_runs[i], seed);
          });
        });
    }
    TBB::instance().wait();

    // Select best partition
    KaHyParPartitioningResult best;
    best.objective = std::numeric_limits<kahypar_hyperedge_weight_t>::max();
    for (size_t i = 0; i < num_ip_calls; ++i) {
      bool improved_metric_within_balance = (results[i].imbalance <= _context.partition.epsilon) &&
                                            (results[i].objective < best.objective);
      bool improved_balance_with_less_equal_metric = (results[i].imbalance < best.imbalance) &&
                                                     (results[i].objective <= best.objective);
      if (improved_metric_within_balance || improved_balance_with_less_equal_metric) {
        best = std::move(results[i]);
      }
    }
    DBG << "Partitioning Result" << V(best.objective) << V(best.imbalance);

    // Apply partition to hypergraph
    for (HypernodeID u = 0; u < best.partition.size(); ++u) {
      ASSERT(u < reverse_mapping.size());
      HypernodeID hn = reverse_mapping[u];
      ASSERT(node_mapping[_hg.originalNodeID(hn)] != kInvalidHypernode);
      PartitionID part_id = best.partition[node_mapping[_hg.originalNodeID(hn)]];
      ASSERT(part_id >= 0 && part_id < _context.partition.k);
      _hg.setNodePart(hn, part_id);
    }
    _hg.updateGlobalPartInfos();
    _hg.initializeNumCutHyperedges();

    ASSERT(metrics::objective(_hg, _context.partition.objective) == best.objective,
      V(metrics::objective(_hg, _context.partition.objective)) << V(best.objective));
    ASSERT(metrics::imbalance(_hg, _context) == best.imbalance);
    kahypar_context_free(context);
  }

 private:
  KaHyParPartitioningResult partitionWithKaHyPar(KaHyParHypergraph& kahypar_hypergraph,
                                                 kahypar_context_t* context,
                                                 const size_t runs,
                                                 const size_t seed) {
    DBG << "Start initial partitioning with" << runs << "initial partitioning runs"
        << "on numa node" << HwTopology::instance().numa_node_of_cpu(sched_getcpu())
        << "on cpu" << sched_getcpu();

    KaHyParPartitioningResult result(kahypar_hypergraph.num_vertices);
    result.objective = std::numeric_limits<kahypar_hyperedge_weight_t>::max();

    // We call KaHyPar exactly one time with runs calls to the initial partitioner of KaHyPar
    // Setup initial partitioning runs
    kahypar::Context initial_partitioning_context(*reinterpret_cast<kahypar::Context*>(context));
    initial_partitioning_context.initial_partitioning.nruns = runs;
    initial_partitioning_context.partition.seed = seed;

    // Call KaHyPar
    kahypar_partition(kahypar_hypergraph, initial_partitioning_context, result);
    result.imbalance = imbalance(kahypar_hypergraph, _context, result);

    DBG << "Finished initial partitioning"
        << "on numa node" << HwTopology::instance().numa_node_of_cpu(sched_getcpu())
        << "on cpu" << sched_getcpu()
        << "with quality" << result.objective
        << "and imbalance" << result.imbalance;
    return result;
  }

  HyperGraph& _hg;
  const Context& _context;
};

template <typename TypeTraits>
PartitionID DirectInitialPartitionerT<TypeTraits>::kInvalidPartition = -1;
template <typename TypeTraits>
HypernodeID DirectInitialPartitionerT<TypeTraits>::kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

using DirectInitialPartitioner = DirectInitialPartitionerT<GlobalTypeTraits>;
}  // namespace mt_kahypar

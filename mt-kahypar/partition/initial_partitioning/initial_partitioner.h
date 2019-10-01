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
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {

template< typename TypeTraits >
class InitialPartitionerT {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;

  static constexpr bool debug = false;
  static constexpr HypernodeID kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

  struct InitialPartitioningResult {

    InitialPartitioningResult() :
      objective(0),
      partition() { }

    explicit InitialPartitioningResult(const kahypar_hypernode_id_t num_vertices) :
      objective(0),
      partition(num_vertices, -1) { }

    kahypar_hyperedge_weight_t objective;
    parallel::scalable_vector<kahypar_partition_id_t> partition;
  };

 public:
  InitialPartitionerT(HyperGraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _context(context) { }

  InitialPartitionerT(const InitialPartitionerT&) = delete;
  InitialPartitionerT(InitialPartitionerT&&) = delete;
  InitialPartitionerT& operator= (const InitialPartitionerT&) = delete;
  InitialPartitionerT& operator= (InitialPartitionerT&&) = delete;

  ~InitialPartitionerT() = default;

  void initialPartition() {
    kahypar_context_t* context = readContext(_context.initial_partitioning.context_file);
    setupContext(*reinterpret_cast<kahypar::Context*>(context));

    // Setup number of runs per thread
    std::vector<size_t> ip_runs;
    if ( _context.initial_partitioning.runs <= _context.shared_memory.num_threads ) {
      ip_runs.assign(_context.initial_partitioning.runs, 1);
    } else {
      size_t runs_per_thread = _context.initial_partitioning.runs / _context.shared_memory.num_threads;
      ip_runs.assign(_context.shared_memory.num_threads, runs_per_thread);
      for ( size_t i = 0; i < _context.initial_partitioning.runs % _context.shared_memory.num_threads; ++i ) {
        ++ip_runs[i];
      }
    }
    ASSERT(ip_runs.size() <= _context.shared_memory.num_threads);
    ASSERT([&] {
      size_t runs = 0;
      for ( const size_t& n : ip_runs ) {
        runs += n;
      }
      return runs == _context.initial_partitioning.runs;
    }(), "Number of runs per thread does not sum up to total number of ip runs");

    // Setup node mapping to a continous range
    HypernodeID num_vertices = 0;
    std::vector<HypernodeID> node_mapping(_hg.initialNumNodes(), kInvalidHypernode);
    for ( const HypernodeID& hn : _hg.nodes() ) {
      ASSERT(_hg.originalNodeID(hn) < _hg.initialNumNodes());
      node_mapping[_hg.originalNodeID(hn)] = num_vertices++;
    }

    // Split calls to initial partitioner evenly across numa nodes
    size_t num_ip_calls = ip_runs.size();
    std::vector<InitialPartitioningResult> results(num_ip_calls);
    size_t used_numa_nodes = TBB::instance().num_used_numa_nodes();
    size_t current_node = 0;
    std::vector<int> used_threads(used_numa_nodes, 0);
    for ( size_t i = 0; i < num_ip_calls; ++i ) {
      while ( used_threads[current_node] >= TBB::instance().number_of_threads_on_numa_node(current_node) ) {
        current_node = (current_node + 1) % used_numa_nodes;
      }

      TBB::instance().numa_task_arena(current_node).execute([&, i] {
        TBB::instance().numa_task_group(current_node).run([&, i] {
          size_t seed = _context.partition.seed + i * _context.initial_partitioning.runs;
          results[i] = partitionWithKaHyPar(context, node_mapping, ip_runs[i], seed);
        });
      });

      ++used_threads[current_node];
      current_node = (current_node + 1) % used_numa_nodes;
    }
    TBB::instance().wait();

    InitialPartitioningResult best;
    best.objective = std::numeric_limits<kahypar_hyperedge_weight_t>::max();
    for ( size_t i = 0; i < num_ip_calls; ++i ) {
      // TODO(heuer): consider imbalance also here
      if ( results[i].objective < best.objective ) {
        best = std::move(results[i]);
      }
    }
    DBG << "Partitioning Result" << V(best.objective);

    // TODO(heuer): Apply best partition to hypergraph

    kahypar_context_free(context);
  }

 private:

  kahypar_context_t* readContext(const std::string& context_file) {
    kahypar_context_t* context = kahypar_context_new();
    kahypar_configure_context_from_file(context, context_file.c_str());
    return context;
  }

  void setupContext(kahypar::Context& context) {
    context.partition.objective = _context.partition.objective;
    context.partition.epsilon = _context.partition.epsilon;
    context.partition.k = _context.partition.k;
    context.partition.seed = _context.partition.seed;
    context.preprocessing.enable_deduplication = true;
    context.preprocessing.enable_min_hash_sparsifier = false;
    context.preprocessing.enable_community_detection = false;
    context.partition.verbose_output = debug;
  }

  InitialPartitioningResult partitionWithKaHyPar(kahypar_context_t* context,
                                                 const std::vector<HypernodeID>& node_mapping,
                                                 const size_t runs,
                                                 const size_t seed) {
    DBG << "Start initial partitioning with" << runs << "initial partitioning runs"
        << "on numa node" << HwTopology::instance().numa_node_of_cpu(sched_getcpu())
        << "on cpu" << sched_getcpu();

    // Setup hypernodes
    kahypar_hypernode_id_t num_vertices = 0;
    parallel::scalable_vector<kahypar_hypernode_weight_t> vertex_weights;
    for ( const HypernodeID& hn : _hg.nodes() ) {
      ASSERT(_hg.originalNodeID(hn) < _hg.initialNumNodes());
      ++num_vertices;
      vertex_weights.emplace_back(_hg.nodeWeight(hn));
    }

    // Setup hyperedges
    kahypar_hyperedge_id_t num_hyperedges = 0;
    parallel::scalable_vector<size_t> hyperedge_indices(1, 0);
    parallel::scalable_vector<kahypar_hyperedge_id_t> hyperedges;
    parallel::scalable_vector<kahypar_hyperedge_weight_t> hyperedge_weights;
    for ( const HyperedgeID& he : _hg.edges() ) {
      ++num_hyperedges;
      for ( const HypernodeID& hn : _hg.pins(he) ) {
        ASSERT(_hg.originalNodeID(hn) < _hg.initialNumNodes());
        ASSERT( node_mapping[_hg.originalNodeID(hn)] != kInvalidHypernode );
        hyperedges.emplace_back(node_mapping[_hg.originalNodeID(hn)]);
      }
      hyperedge_indices.emplace_back(hyperedges.size());
      hyperedge_weights.emplace_back(_hg.edgeWeight(he));
    }

    InitialPartitioningResult best(num_vertices);
    best.objective = std::numeric_limits<kahypar_hyperedge_weight_t>::max();
    if ( _context.initial_partitioning.call_kahypar_multiple_times ) {
      // We call KaHyPar exactly runs times with one call to the initial partitioner of KaHyPar
      for ( size_t i = 0; i < runs; ++i ) {
        // Setup initial partitioning runs
        kahypar::Context initial_partitioning_context(*reinterpret_cast<kahypar::Context*>(context));
        initial_partitioning_context.initial_partitioning.nruns = 1;
        initial_partitioning_context.partition.seed = seed + i;

        // Call KaHyPar
        InitialPartitioningResult result(num_vertices);
        kahypar_partition(num_vertices, num_hyperedges, _context.partition.epsilon,
          _context.partition.k, vertex_weights.data(), hyperedge_weights.data(),
          hyperedge_indices.data(), hyperedges.data(), &result.objective,
          reinterpret_cast<kahypar_context_t*>(&initial_partitioning_context),
          result.partition.data());

        // TODO(heuer): consider imbalance also here
        if ( result.objective < best.objective ) {
          best = std::move(result);
        }
      }
    } else {
      // We call KaHyPar exactly one time with runs calls to the initial partitioner of KaHyPar
      // Setup initial partitioning runs
      kahypar::Context initial_partitioning_context(*reinterpret_cast<kahypar::Context*>(context));
      initial_partitioning_context.initial_partitioning.nruns = runs;
      initial_partitioning_context.partition.seed = seed;

      // Call KaHyPar
      kahypar_partition(num_vertices, num_hyperedges, _context.partition.epsilon,
        _context.partition.k, vertex_weights.data(), hyperedge_weights.data(),
        hyperedge_indices.data(), hyperedges.data(), &best.objective,
        reinterpret_cast<kahypar_context_t*>(&initial_partitioning_context),
        best.partition.data());
    }

    DBG << "Finish initial partitioning"
        << "on numa node" << HwTopology::instance().numa_node_of_cpu(sched_getcpu())
        << "on cpu" << sched_getcpu()
        << "with quality" << best.objective;
    return best;
  }

  HyperGraph& _hg;
  const Context& _context;
};

using InitialPartitioner = InitialPartitionerT<GlobalTypeTraits>;

}  // namespace mt_kahypar

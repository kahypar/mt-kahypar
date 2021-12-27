/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "libmtkahypar.h"

#include "tbb/parallel_for.h"
#include "tbb/parallel_invoke.h"

#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/partition/partitioner.h"
#include "mt-kahypar/partition/registries/register_memory_pool.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/utils/randomize.h"

namespace {
  template<typename T>
  using vec = mt_kahypar::parallel::scalable_vector<T>;
}

mt_kahypar_context_t* mt_kahypar_context_new() {
  return reinterpret_cast<mt_kahypar_context_t*>(new mt_kahypar::Context());
}

void mt_kahypar_context_free(mt_kahypar_context_t* kahypar_context) {
  if (kahypar_context == nullptr) {
    return;
  }
  delete reinterpret_cast<mt_kahypar::Context*>(kahypar_context);
}

void mt_kahypar_configure_context_from_file(mt_kahypar_context_t* kahypar_context,
                                            const char* ini_file_name) {
  mt_kahypar::parseIniToContext(*reinterpret_cast<mt_kahypar::Context*>(kahypar_context),
                                ini_file_name);
}


void mt_kahypar_initialize_thread_pool(const size_t num_threads,
                                       const bool interleaved_allocations) {
  size_t P = num_threads;
  size_t num_available_cpus = mt_kahypar::HardwareTopology::instance().num_cpus();
  if ( num_available_cpus < num_threads ) {
    WARNING("There are currently only" << num_available_cpus << "cpus available."
      << "Setting number of threads from" << num_threads
      << "to" << num_available_cpus);
    P = num_available_cpus;
  }

  // Initialize TBB task arenas on numa nodes
  mt_kahypar::TBBInitializer::instance(P);

  if ( interleaved_allocations ) {
    // We set the membind policy to interleaved allocations in order to
    // distribute allocations evenly across NUMA nodes
    hwloc_cpuset_t cpuset = mt_kahypar::TBBInitializer::instance().used_cpuset();
    mt_kahypar::parallel::HardwareTopology<>::instance().activate_interleaved_membind_policy(cpuset);
    hwloc_bitmap_free(cpuset);
  }
}


void mt_kahypar_read_hypergraph_from_file(const char* file_name,
                                          mt_kahypar_hypernode_id_t* num_vertices,
                                          mt_kahypar_hyperedge_id_t* num_hyperedges,
                                          size_t** hyperedge_indices,
                                          mt_kahypar_hyperedge_id_t** hyperedges,
                                          mt_kahypar_hyperedge_weight_t** hyperedge_weights,
                                          mt_kahypar_hypernode_weight_t** vertex_weights) {
  mt_kahypar::HypernodeID num_nodes = 0;
  mt_kahypar::HyperedgeID num_edges = 0;
  mt_kahypar::HyperedgeID num_removed_hyperedges = 0;
  vec<vec<mt_kahypar::HypernodeID>> edge_vector;
  vec<mt_kahypar::HypernodeWeight> hypernodes_weight;
  vec<mt_kahypar::HyperedgeWeight> net_weight;
  // TODO(maas): support for other file formats
  mt_kahypar::io::readHypergraphFile(file_name, num_edges, num_nodes,
    num_removed_hyperedges, edge_vector, net_weight, hypernodes_weight);
  ASSERT(num_edges == edge_vector.size());
  ASSERT(num_edges == net_weight.size());
  ASSERT(num_nodes == hypernodes_weight.size());

  *num_vertices = num_nodes;
  *num_hyperedges = num_edges;
  tbb::parallel_invoke([&] {
    *hyperedge_indices = new size_t[num_edges + 1];

    // Compute index vector with parallel prefix sum
    vec<size_t> indices(num_edges + 1, 0);
    tbb::parallel_for(0UL, static_cast<size_t>(num_edges), [&](const size_t i) {
      indices[i + 1] = edge_vector[i].size();
    });
    mt_kahypar::parallel::TBBPrefixSum<size_t, vec> index_prefix_sum(indices);
    tbb::parallel_scan(tbb::blocked_range<size_t>(
      0UL, UI64(num_edges + 1)), index_prefix_sum);

    // Setup Incidence Array
    tbb::parallel_invoke([&] {
      *hyperedges = new mt_kahypar_hyperedge_id_t[index_prefix_sum.total_sum()];
      tbb::parallel_for(0UL, static_cast<size_t>(num_edges), [&](const size_t i) {
        for ( size_t j = 0; j < edge_vector[i].size(); ++j ) {
          ASSERT(indices[i] + j < indices[i + 1]);
          *(*hyperedges + (indices[i] + j)) = edge_vector[i][j];
        }
      });
    }, [&] {
      memcpy(*hyperedge_indices, indices.data(), sizeof(size_t) * (num_edges + 1));
    });
  }, [&] {
    *hyperedge_weights = new mt_kahypar_hyperedge_weight_t[num_edges];
    memcpy(*hyperedge_weights, net_weight.data(),
      sizeof(mt_kahypar_hyperedge_weight_t) * num_edges);
  }, [&] {
    *vertex_weights = new mt_kahypar_hypernode_weight_t[num_nodes];
    memcpy(*vertex_weights, hypernodes_weight.data(),
      sizeof(mt_kahypar_hypernode_weight_t) * num_nodes);
  });
}

void mt_kahypar_partition(const mt_kahypar_hypernode_id_t num_vertices,
                          const mt_kahypar_hyperedge_id_t num_hyperedges,
                          const double epsilon,
                          const mt_kahypar_partition_id_t num_blocks,
                          const int seed,
                          const mt_kahypar_hypernode_weight_t* vertex_weights,
                          const mt_kahypar_hyperedge_weight_t* hyperedge_weights,
                          const size_t* hyperedge_indices,
                          const mt_kahypar_hyperedge_id_t* hyperedges,
                          mt_kahypar_hyperedge_weight_t* objective,
                          mt_kahypar_context_t* kahypar_context,
                          mt_kahypar_partition_id_t* partition,
                          const bool verbose) {
  mt_kahypar::Context context = *reinterpret_cast<mt_kahypar::Context*>(kahypar_context);
  context.partition.k = num_blocks;
  context.partition.epsilon = epsilon;
  context.partition.mode = mt_kahypar::Mode::direct;
  context.partition.objective = kahypar::Objective::km1;
  context.partition.seed = seed;
  context.partition.verbose_output = verbose;
  context.partition.write_partition_file = false;

  mt_kahypar::utils::Randomize::instance().setSeed(context.partition.seed);

  // TODO(heuer): change internal hypergraph construction format
  // Transform adjacence array into adjacence list
  vec<vec<mt_kahypar::HypernodeID>> edge_vector(num_hyperedges);
  tbb::parallel_for(0UL, num_hyperedges, [&](const mt_kahypar::HyperedgeID& he) {
    const size_t num_pins = hyperedge_indices[he + 1] - hyperedge_indices[he];
    edge_vector[he].resize(num_pins);
    for ( size_t i = 0; i < num_pins; ++i ) {
      edge_vector[he][i] = hyperedges[hyperedge_indices[he] + i];
    }
  });

  // Contruct Hypergraph
  mt_kahypar::Hypergraph hypergraph = mt_kahypar::HypergraphFactory::construct(
    num_vertices, num_hyperedges, edge_vector, hyperedge_weights, vertex_weights);

  // Partition Hypergraph
  mt_kahypar::PartitionedHypergraph partitioned_hypergraph =
    mt_kahypar::partition(hypergraph, context);

  // Store partition
  *objective = mt_kahypar::metrics::objective(partitioned_hypergraph, context.partition.objective);
  ASSERT(partition != nullptr);
  partitioned_hypergraph.doParallelForAllNodes([&](const mt_kahypar::HypernodeID& hn) {
    partition[hn] = partitioned_hypergraph.partID(hn);
  });
}
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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/pin_count_in_part.h"
#include "mt-kahypar/datastructures/connectivity_set.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/memory_tree.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

static void register_memory_pool(const Hypergraph& hypergraph,
                                 const Context& context) {

  // ########## Preprocessing Memory ##########

  if ( context.preprocessing.use_community_detection ) {
    const HypernodeID num_hypernodes = hypergraph.initialNumNodes();
    const HyperedgeID num_hyperedges = hypergraph.initialNumEdges();
    const HypernodeID num_pins = hypergraph.initialNumPins();
    const bool is_graph = hypergraph.maxEdgeSize() == 2;
    const size_t num_nodes = num_hypernodes + (is_graph ? 0 : num_hyperedges);
    const size_t num_edges = is_graph ? num_pins : (2UL * num_pins);

    parallel::MemoryPool::instance().register_memory_group("Preprocessing", 1);
    parallel::MemoryPool::instance().register_memory_chunk("Preprocessing", "indices", num_nodes + 1, sizeof(size_t));
    parallel::MemoryPool::instance().register_memory_chunk("Preprocessing", "arcs", num_edges, sizeof(Arc));
    parallel::MemoryPool::instance().register_memory_chunk("Preprocessing", "node_volumes", num_nodes, sizeof(ArcWeight));
    parallel::MemoryPool::instance().register_memory_chunk("Preprocessing", "tmp_indices",
      num_nodes + 1, sizeof(parallel::IntegralAtomicWrapper<size_t>));
    parallel::MemoryPool::instance().register_memory_chunk("Preprocessing", "tmp_pos",
      num_nodes, sizeof(parallel::IntegralAtomicWrapper<size_t>));
    parallel::MemoryPool::instance().register_memory_chunk("Preprocessing", "tmp_arcs", num_edges, sizeof(Arc));
    parallel::MemoryPool::instance().register_memory_chunk("Preprocessing", "valid_arcs", num_edges, sizeof(size_t));
    parallel::MemoryPool::instance().register_memory_chunk("Preprocessing", "tmp_node_volumes",
      num_nodes, sizeof(parallel::AtomicWrapper<ArcWeight>));
  }

  // ########## Coarsening Memory ##########

  parallel::MemoryPool::instance().register_memory_group("Coarsening", 2);
  const bool is_numa_aware = Hypergraph::is_numa_aware;
  for ( size_t node = 0; node < hypergraph.numNumaHypergraphs(); ++node ) {
    const HypernodeID num_hypernodes = hypergraph.initialNumNodes(node);
    const HyperedgeID num_hyperedges = hypergraph.initialNumEdges(node);
    const HypernodeID num_pins = hypergraph.initialNumPins(node);
    parallel::MemoryPool::instance().register_memory_chunk(
      "Coarsening", "mapping_" + std::to_string(node), num_hypernodes, sizeof(size_t), node);
    parallel::MemoryPool::instance().register_memory_chunk(
      "Coarsening", "tmp_hypernodes_" + std::to_string(node), num_hypernodes, Hypergraph::SIZE_OF_HYPERNODE, node);
    if ( !is_numa_aware ) {
      parallel::MemoryPool::instance().register_memory_chunk(
        "Coarsening", "tmp_incident_nets_" + std::to_string(node), num_pins, sizeof(HyperedgeID), node);
    }
    parallel::MemoryPool::instance().register_memory_chunk(
      "Coarsening", "tmp_num_incident_nets_" + std::to_string(node),
      num_hypernodes, sizeof(parallel::IntegralAtomicWrapper<size_t>), node);
    parallel::MemoryPool::instance().register_memory_chunk(
      "Coarsening", "hn_weights_" + std::to_string(node),
      num_hypernodes, sizeof(parallel::IntegralAtomicWrapper<HypernodeWeight>), node);
    parallel::MemoryPool::instance().register_memory_chunk(
      "Coarsening", "tmp_hyperedges_" + std::to_string(node), num_hyperedges, Hypergraph::SIZE_OF_HYPEREDGE, node);
    parallel::MemoryPool::instance().register_memory_chunk(
      "Coarsening", "tmp_incidence_array_" + std::to_string(node), num_pins, sizeof(HypernodeID), node);
    parallel::MemoryPool::instance().register_memory_chunk(
      "Coarsening", "he_sizes_" + std::to_string(node), num_hyperedges, sizeof(size_t), node);
    parallel::MemoryPool::instance().register_memory_chunk(
      "Coarsening", "valid_hyperedges_" + std::to_string(node), num_hyperedges, sizeof(size_t), node);
  }

  // ########## Refinement Memory ##########

  parallel::MemoryPool::instance().register_memory_group("Refinement", 3);

  for ( size_t node = 0; node < hypergraph.numNumaHypergraphs(); ++node ) {
    const HypernodeID num_hypernodes = hypergraph.initialNumNodes(node);
    const HyperedgeID num_hyperedges = hypergraph.initialNumEdges(node);
    const HypernodeID max_he_size = hypergraph.maxEdgeSize();
    parallel::MemoryPool::instance().register_memory_chunk(
      "Refinement", "vertex_part_info_" + std::to_string(node),
      num_hypernodes, PartitionedHypergraph<>::SIZE_OF_VERTEX_PART_INFO, node);
    parallel::MemoryPool::instance().register_memory_chunk(
      "Refinement", "pin_count_in_part_" + std::to_string(node),
      ds::PinCountInPart::num_elements(num_hyperedges, context.partition.k, max_he_size),
      sizeof(ds::PinCountInPart::Value), node);
    parallel::MemoryPool::instance().register_memory_chunk(
      "Refinement", "connectivity_set_" + std::to_string(node),
      ds::ConnectivitySets::num_elements(num_hyperedges, context.partition.k),
      sizeof(ds::ConnectivitySets::UnsafeBlock), node);
    parallel::MemoryPool::instance().register_memory_chunk(
      "Refinement", "pin_count_update_ownership_" + std::to_string(node),
      num_hyperedges, sizeof(parallel::IntegralAtomicWrapper<bool>), node);
  }

  // Allocate Memory
  utils::Timer::instance().start_timer("memory_pool_allocation", "Memory Pool Allocation");
  parallel::MemoryPool::instance().allocate_memory_chunks<TBBNumaArena>();
  utils::Timer::instance().stop_timer("memory_pool_allocation");
}

} // namespace mt_kahypar
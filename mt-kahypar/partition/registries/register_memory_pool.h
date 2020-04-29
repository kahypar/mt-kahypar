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
  const HypernodeID num_hypernodes = hypergraph.initialNumNodes();
  const HyperedgeID num_hyperedges = hypergraph.initialNumEdges();
  const HypernodeID num_pins = hypergraph.initialNumPins();
  parallel::MemoryPool::instance().register_memory_chunk("Coarsening", "mapping", num_hypernodes, sizeof(size_t));
  parallel::MemoryPool::instance().register_memory_chunk("Coarsening", "tmp_hypernodes", num_hypernodes, Hypergraph::SIZE_OF_HYPERNODE);
  parallel::MemoryPool::instance().register_memory_chunk("Coarsening", "tmp_incident_nets", num_pins, sizeof(HyperedgeID));
  parallel::MemoryPool::instance().register_memory_chunk("Coarsening", "tmp_num_incident_nets",
    num_hypernodes, sizeof(parallel::IntegralAtomicWrapper<size_t>));
  parallel::MemoryPool::instance().register_memory_chunk("Coarsening", "hn_weights",
    num_hypernodes, sizeof(parallel::IntegralAtomicWrapper<HypernodeWeight>));
  parallel::MemoryPool::instance().register_memory_chunk("Coarsening", "tmp_hyperedges", num_hyperedges, Hypergraph::SIZE_OF_HYPEREDGE);
  parallel::MemoryPool::instance().register_memory_chunk("Coarsening", "tmp_incidence_array", num_pins, sizeof(HypernodeID));
  parallel::MemoryPool::instance().register_memory_chunk("Coarsening", "he_sizes", num_hyperedges, sizeof(size_t));
  parallel::MemoryPool::instance().register_memory_chunk("Coarsening", "valid_hyperedges", num_hyperedges, sizeof(size_t));

  // ########## Refinement Memory ##########

  parallel::MemoryPool::instance().register_memory_group("Refinement", 3);
  const HypernodeID max_he_size = hypergraph.maxEdgeSize();
  parallel::MemoryPool::instance().register_memory_chunk("Refinement", "part_ids", num_hypernodes, sizeof(PartitionID));
  parallel::MemoryPool::instance().register_memory_chunk("Refinement", "pin_count_in_part",
    ds::PinCountInPart::num_elements(num_hyperedges, context.partition.k, max_he_size),
    sizeof(ds::PinCountInPart::Value));
  parallel::MemoryPool::instance().register_memory_chunk("Refinement", "connectivity_set",
    ds::ConnectivitySets::num_elements(num_hyperedges, context.partition.k),
    sizeof(ds::ConnectivitySets::UnsafeBlock));
  parallel::MemoryPool::instance().register_memory_chunk("Refinement", "move_to_penalty",
    num_hypernodes * context.partition.k, sizeof(CAtomic<HyperedgeWeight>));
  parallel::MemoryPool::instance().register_memory_chunk("Refinement", "move_from_penalty",
    num_hypernodes, sizeof(CAtomic<HyperedgeWeight>));
  parallel::MemoryPool::instance().register_memory_chunk("Refinement", "pin_count_update_ownership",
    num_hyperedges, sizeof(parallel::IntegralAtomicWrapper<bool>));

  // Allocate Memory
  utils::Timer::instance().start_timer("memory_pool_allocation", "Memory Pool Allocation");
  parallel::MemoryPool::instance().allocate_memory_chunks<TBBNumaArena>();
  utils::Timer::instance().stop_timer("memory_pool_allocation");
}

} // namespace mt_kahypar
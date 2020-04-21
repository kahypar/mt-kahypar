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
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/memory_tree.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

static void register_memory_pool(const Hypergraph& hypergraph) {

  const HypernodeID num_hypernodes = hypergraph.initialNumNodes();
  const HyperedgeID num_hyperedges = hypergraph.initialNumEdges();
  const HypernodeID num_pins = hypergraph.initialNumPins();

  // ########## Community Detection Memory ##########
  const bool is_graph = hypergraph.maxEdgeSize() == 2;
  const size_t num_nodes = num_hypernodes + (is_graph ? 0 : num_hyperedges);
  const size_t num_edges = is_graph ? num_pins : (2UL * num_pins);

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

  // Allocate Memory
  utils::Timer::instance().start_timer("memory_pool_allocation", "Memory Pool Allocation");
  parallel::MemoryPool::instance().allocate_memory_chunks();
  utils::Timer::instance().stop_timer("memory_pool_allocation");
}

} // namespace mt_kahypar
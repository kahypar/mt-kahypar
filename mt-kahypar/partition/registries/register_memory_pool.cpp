/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include "register_memory_pool.h"

#include "mt-kahypar/datastructures/pin_count_in_part.h"
#include "mt-kahypar/datastructures/connectivity_set.h"
#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/memory_tree.h"
#include "mt-kahypar/utils/utilities.h"

namespace mt_kahypar {

  void register_memory_pool(const Hypergraph& hypergraph,
                            const Context& context) {

    if (context.partition.mode == Mode::direct) {

      // ########## Preprocessing Memory ##########

      const HypernodeID num_hypernodes = hypergraph.initialNumNodes();
      const HyperedgeID num_hyperedges = hypergraph.initialNumEdges();
      const HypernodeID num_pins = hypergraph.initialNumPins();

      auto& pool = parallel::MemoryPool::instance();

      if ( context.preprocessing.use_community_detection ) {
        const bool is_graph = hypergraph.maxEdgeSize() == 2;
        const size_t num_star_expansion_nodes = num_hypernodes + (is_graph ? 0 : num_hyperedges);
        const size_t num_star_expansion_edges = is_graph ? num_pins : (2UL * num_pins);

        pool.register_memory_group("Preprocessing", 1);
        pool.register_memory_chunk("Preprocessing", "indices", num_star_expansion_nodes + 1, sizeof(size_t));
        pool.register_memory_chunk("Preprocessing", "arcs", num_star_expansion_edges, sizeof(Arc));
        pool.register_memory_chunk("Preprocessing", "node_volumes", num_star_expansion_nodes, sizeof(ArcWeight));

        if ( !context.preprocessing.community_detection.low_memory_contraction ) {
          pool.register_memory_chunk("Preprocessing", "tmp_indices",
                                    num_star_expansion_nodes + 1, sizeof(parallel::IntegralAtomicWrapper<size_t>));
          pool.register_memory_chunk("Preprocessing", "tmp_pos",
                                    num_star_expansion_nodes, sizeof(parallel::IntegralAtomicWrapper<size_t>));
          pool.register_memory_chunk("Preprocessing", "tmp_arcs", num_star_expansion_edges, sizeof(Arc));
          pool.register_memory_chunk("Preprocessing", "valid_arcs", num_star_expansion_edges, sizeof(size_t));
          pool.register_memory_chunk("Preprocessing", "tmp_node_volumes",
                                    num_star_expansion_nodes, sizeof(parallel::AtomicWrapper<ArcWeight>));
        }
      }

      // ########## Coarsening Memory ##########

      pool.register_memory_group("Coarsening", 2);
      if ( context.partition.paradigm == Paradigm::multilevel ) {
        if (Hypergraph::is_graph) {
          pool.register_memory_chunk("Coarsening", "mapping", num_hypernodes, sizeof(HypernodeID));
          pool.register_memory_chunk("Coarsening", "tmp_nodes", num_hypernodes, Hypergraph::SIZE_OF_HYPERNODE);
          pool.register_memory_chunk("Coarsening", "node_sizes", num_hypernodes, sizeof(HyperedgeID));
          pool.register_memory_chunk("Coarsening", "tmp_num_incident_edges",
                                     num_hypernodes, sizeof(parallel::IntegralAtomicWrapper<HyperedgeID>));
          pool.register_memory_chunk("Coarsening", "node_weights",
                                     num_hypernodes, sizeof(parallel::IntegralAtomicWrapper<HypernodeWeight>));
          pool.register_memory_chunk("Coarsening", "tmp_edges", num_hyperedges, Hypergraph::SIZE_OF_HYPEREDGE);
          pool.register_memory_chunk("Coarsening", "edge_id_mapping", num_hyperedges / 2, sizeof(HyperedgeID));
        } else {
          pool.register_memory_chunk("Coarsening", "mapping", num_hypernodes, sizeof(size_t));
          pool.register_memory_chunk("Coarsening", "tmp_hypernodes", num_hypernodes, Hypergraph::SIZE_OF_HYPERNODE);
          pool.register_memory_chunk("Coarsening", "tmp_incident_nets", num_pins, sizeof(HyperedgeID));
          pool.register_memory_chunk("Coarsening", "tmp_num_incident_nets",
                                    num_hypernodes, sizeof(parallel::IntegralAtomicWrapper<size_t>));
          pool.register_memory_chunk("Coarsening", "hn_weights",
                                    num_hypernodes, sizeof(parallel::IntegralAtomicWrapper<HypernodeWeight>));
          pool.register_memory_chunk("Coarsening", "tmp_hyperedges", num_hyperedges, Hypergraph::SIZE_OF_HYPEREDGE);
          pool.register_memory_chunk("Coarsening", "tmp_incidence_array", num_pins, sizeof(HypernodeID));
          pool.register_memory_chunk("Coarsening", "he_sizes", num_hyperedges, sizeof(size_t));
          pool.register_memory_chunk("Coarsening", "valid_hyperedges", num_hyperedges, sizeof(size_t));
        }
      }

      // ########## Refinement Memory ##########

      pool.register_memory_group("Refinement", 3);
      pool.register_memory_chunk("Refinement", "part_ids", num_hypernodes, sizeof(PartitionID));

      if (Hypergraph::is_graph) {
        #ifdef USE_GRAPH_PARTITIONER // SIZE_OF_EDGE_LOCK is only available in the graph data structure
          pool.register_memory_chunk("Refinement", "edge_locks", num_hyperedges, PartitionedHypergraph::SIZE_OF_EDGE_LOCK);
        #endif
        if ( context.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
          pool.register_memory_chunk("Refinement", "incident_weight_in_part",
                                    static_cast<size_t>(num_hypernodes) * ( context.partition.k + 1 ),
                                    sizeof(CAtomic<HyperedgeWeight>));
        }
      } else {
        const HypernodeID max_he_size = hypergraph.maxEdgeSize();
        pool.register_memory_chunk("Refinement", "pin_count_in_part",
                                  ds::PinCountInPart::num_elements(num_hyperedges, context.partition.k, max_he_size),
                                  sizeof(ds::PinCountInPart::Value));
        pool.register_memory_chunk("Refinement", "connectivity_set",
                                  ds::ConnectivitySets::num_elements(num_hyperedges, context.partition.k),
                                  sizeof(ds::ConnectivitySets::UnsafeBlock));
        if ( context.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
          pool.register_memory_chunk("Refinement", "gain_cache",
                                    static_cast<size_t>(num_hypernodes) * ( context.partition.k + 1 ),
                                    sizeof(CAtomic<HyperedgeWeight>));
        }
        pool.register_memory_chunk("Refinement", "pin_count_update_ownership",
                                  num_hyperedges, sizeof(SpinLock));
      }

      // Allocate Memory
      utils::Timer& timer = utils::Utilities::instance().getTimer(context.utility_id);
      timer.start_timer("memory_pool_allocation", "Memory Pool Allocation");
      pool.allocate_memory_chunks();
      timer.stop_timer("memory_pool_allocation");
    }
  }


} // namespace mt_kahypar

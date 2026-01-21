/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
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

#pragma once

#include "mt-kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/initial_partitioning_data_container.h"
#include "mt-kahypar/parallel/stl/scalable_queue.h"
#include "mt-kahypar/partition/constraints.h"
#include "mt-kahypar/partition/k_colouring.h"

namespace mt_kahypar {

template<typename TypeTraits>
class BFSInitialPartitioner : public IInitialPartitioner {
  using Queue = parallel::scalable_queue<HypernodeID>;

  static constexpr bool debug = false;

  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;

  struct transformed_colouring {
    vec<Colour> node_colours;
    vec<HypernodeID> nodes_per_colour;
    Colour current_max_colour;

    template<typename PartitionedHypergraph>
    static transformed_colouring get_colouring(const PartitionedHypergraph& hypergraph, const Context& context) {
      graph_colouring constraint_colouring = KColouring<PartitionedHypergraph>(context).colour(hypergraph);
      transformed_colouring transformed;
      transformed.current_max_colour = constraint_colouring.used_colours - 1;
      transformed.nodes_per_colour = vec<HypernodeID> (constraint_colouring.used_colours, 0);
      transformed.node_colours = vec<Colour> (hypergraph.initialNumNodes(), kInvalidColour);
      for (const HypernodeID hn : hypergraph.nodes()) {
        HypernodeID node;
        if (hypergraph.fixedVertexSupport().getConstraintIdFromHypergraphId(hn, node)) {
          Colour colour = constraint_colouring.node_colours[node];
          transformed.nodes_per_colour[colour]++;
          transformed.node_colours[hn] = colour;
        }
      }
      return transformed;
    }

    bool isNodeAllowed(const HypernodeID hn) const {
      Colour colour = node_colours[hn];
      if (colour == kInvalidColour) return true;
      if (colour == current_max_colour) return true;
      return false;
    }

    void setNodePart(const HypernodeID hn, const PartitionID) {
      Colour colour = node_colours[hn];
      if (colour == kInvalidColour) return;
      if (colour != current_max_colour) {
        LOG << "Used not the biggest colour!"<< colour << current_max_colour;
      }
      nodes_per_colour[colour]--;
      if (nodes_per_colour[colour] <= 0 && colour == current_max_colour) {
        current_max_colour--;
      }
    }

  };

 public:
  BFSInitialPartitioner(const InitialPartitioningAlgorithm,
                        ip_data_container_t* ip_data,
                        const Context& context,
                        const int seed, const int tag) :
    _ip_data(ip::to_reference<TypeTraits>(ip_data)),
    _context(context),
    _rng(seed),
    _tag(tag) { }

 private:
  void partitionImpl() final;

  bool fitsIntoBlock(PartitionedHypergraph& hypergraph,
                     const HypernodeID hn,
                     const PartitionID block) const {
    ASSERT(block != kInvalidPartition && block < _context.partition.k);
    return hypergraph.partWeight(block) + hypergraph.nodeWeight(hn) <=
      _context.partition.perfect_balance_part_weights[block];
  }

  bool constraintsAllowBlock(PartitionedHypergraph& hypergraph,
                            const HypernodeID hn,
                            const PartitionID block,
                            const transformed_colouring& colouring) const {
    if (hypergraph.hasNegativeConstraints()) {
      if (colouring.isNodeAllowed(hn)) {
        return constraints::isNodeAllowedInPartition(hypergraph, _context, hn, block);
      }
      return false;
    }
    return true;
  }

  bool isConstraintNode(PartitionedHypergraph& hypergraph,
                        const HypernodeID hn) {
    HypernodeID node;
    return hypergraph.fixedVertexSupport().getConstraintIdFromHypergraphId(hn, node);
  }

  // ! Pushes all adjacent hypernodes (not visited before) of hypernode hn
  // ! into the BFS queue of the corresponding block.
  inline void pushIncidentHypernodesIntoQueue(const PartitionedHypergraph& hypergraph,
                                              const Context& context,
                                              Queue& queue,
                                              kahypar::ds::FastResetFlagArray<>& hypernodes_in_queue,
                                              kahypar::ds::FastResetFlagArray<>& hyperedges_in_queue,
                                              const HypernodeID hn,
                                              const PartitionID block);

  inline void markHypernodeAsInQueue(const PartitionedHypergraph& hypergraph,
                                     kahypar::ds::FastResetFlagArray<>& hypernodes_in_queue,
                                     const HypernodeID hn,
                                     const PartitionID block) {
    hypernodes_in_queue.set(block * hypergraph.initialNumNodes() + hn, true);
  }

  inline void markHyperedgeAsInQueue(const PartitionedHypergraph& hypergraph,
                                     kahypar::ds::FastResetFlagArray<>& hyperedges_in_queue,
                                     const HyperedgeID he,
                                     const PartitionID block) {
    hyperedges_in_queue.set(block * hypergraph.initialNumEdges() + he, true);
  }

  InitialPartitioningDataContainer<TypeTraits>& _ip_data;
  const Context& _context;
  std::mt19937 _rng;
  const int _tag;
};


} // namespace mt_kahypar

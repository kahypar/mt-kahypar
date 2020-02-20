/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include <vector>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/clustering.h"
#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {
class HypergraphSparsifier {
 public:

  HypergraphSparsifier() :
    _removed_hes(),
    _removed_hns() { }

  HypergraphSparsifier(const HypergraphSparsifier&) = delete;
  HypergraphSparsifier & operator= (const HypergraphSparsifier &) = delete;

  HypergraphSparsifier(HypergraphSparsifier&&) = delete;
  HypergraphSparsifier & operator= (HypergraphSparsifier &&) = delete;

  HyperedgeID removeSingleNodeHyperedges(Hypergraph& hypergraph) {
    HyperedgeID num_removed_single_node_hes = 0;
    for (const HyperedgeID& he : hypergraph.edges()) {
      if (hypergraph.edgeSize(he) == 1) {
        ++num_removed_single_node_hes;
        hypergraph.removeEdge(he);
        _removed_hes.push_back(he);
      }
    }
    return num_removed_single_node_hes;
  }

  HypernodeID removeDegreeZeroHypernodes(Hypergraph& hypergraph) {
    HypernodeID num_removed_degree_zero_hypernodes = 0;
    for (const HypernodeID& hn : hypergraph.nodes()) {
      // Currently, we only remove zero-degree hypernodes with weight one
      // If, we would remove hypernodes with weight greater than one than a
      // a feasible partition would be not always possible, if we restore
      // them afterwards.
      if ( hypergraph.nodeDegree(hn) == 0 && hypergraph.nodeWeight(hn) == 1) {
        ++num_removed_degree_zero_hypernodes;
        hypergraph.removeHypernode(hn);
        _removed_hns.push_back(hn);
      }
    }
    return num_removed_degree_zero_hypernodes;
  }

  void assignAllDegreeZeroHypernodesToSameCommunity(Hypergraph& hypergraph, ds::Clustering& clustering) {
    ASSERT(hypergraph.initialNumNodes() == clustering.size());
    PartitionID community_id = -1;
    for ( const HypernodeID& hn : hypergraph.nodes() ) {
      if ( hypergraph.nodeDegree(hn) == 0 ) {
        if ( community_id >= 0 ) {
          clustering[hypergraph.originalNodeID(hn)] = community_id;
        } else {
          community_id = clustering[hypergraph.originalNodeID(hn)];
        }
      }
    }
  }

  void restoreSingleNodeHyperedges(PartitionedHypergraph<>& hypergraph) {
    for (const HyperedgeID& he : _removed_hes) {
      hypergraph.restoreSinglePinHyperedge(he);
    }
  }

  void restoreDegreeZeroHypernodes(PartitionedHypergraph<>& hypergraph, const Context& context) {
    std::vector<PartitionID> valid_blocks(context.partition.k, 0);
    std::iota(valid_blocks.begin(), valid_blocks.end(), 0);
    size_t current_block_idx = 0;
    for ( const HypernodeID& hn : _removed_hns ) {
      hypergraph.enableHypernode(hn);
      ASSERT(hypergraph.nodeWeight(hn) == 1);

      // Search for a valid block to which we assign the restored hypernode
      // Note, a valid block must exist because we only remove weight one hypernodes
      PartitionID current_block = valid_blocks[current_block_idx];
      while ( hypergraph.partWeight(current_block) >= context.partition.max_part_weights[current_block] ) {
        std::swap(valid_blocks[current_block_idx], valid_blocks.back());
        valid_blocks.pop_back();
        ASSERT(!valid_blocks.empty());
        current_block_idx = current_block_idx % valid_blocks.size();
        current_block = valid_blocks[current_block_idx];
      }

      hypergraph.setNodePart(hn, current_block);
      current_block_idx = (current_block_idx + 1) % valid_blocks.size();
    }
  }

 private:
  std::vector<HyperedgeID> _removed_hes;
  std::vector<HypernodeID> _removed_hns;
};
}  // namespace mt_kahypar

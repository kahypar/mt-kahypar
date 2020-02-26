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
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
class HypergraphSparsifier {
 public:

  HypergraphSparsifier() :
    _removed_hes(),
    _removed_hns(),
    _representative_hns() { }

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

  HypernodeID contractDegreeZeroHypernodes(Hypergraph& hypergraph, const Context& context) {
    _representative_hns.assign(hypergraph.initialNumNodes(), kInvalidHypernode);
    HypernodeID num_removed_degree_zero_hypernodes = 0;
    HypernodeID last_representative = kInvalidHypernode;
    HypernodeWeight last_representative_weight = 0;
    for (const HypernodeID& hn : hypergraph.nodes()) {
      if ( hypergraph.nodeDegree(hn) == 0 ) {
        bool was_removed = false;
        if ( last_representative != kInvalidHypernode ) {
          const HypernodeWeight weight = hypergraph.nodeWeight(hn);
          if ( last_representative_weight + weight <= context.coarsening.max_allowed_node_weight ) {
            ++num_removed_degree_zero_hypernodes;
            hypergraph.removeHypernode(hn);
            hypergraph.setNodeWeight(last_representative, last_representative_weight + weight);
            last_representative_weight += weight;
            _removed_hns.push_back(hn);
            _representative_hns[hypergraph.originalNodeID(hn)] = last_representative;
            was_removed = true;
          }
        }

        if ( !was_removed ) {
          last_representative = hn;
          last_representative_weight = hypergraph.nodeWeight(hn);
        }
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

  void restoreDegreeZeroHypernodes(PartitionedHypergraph<>& hypergraph) {
    for ( const HypernodeID& hn : _removed_hns ) {
      const HypernodeID original_id = hypergraph.originalNodeID(hn);
      ASSERT(original_id < _representative_hns.size());
      const HypernodeID representative = _representative_hns[original_id];
      ASSERT(representative != kInvalidHypernode);
      hypergraph.enableHypernode(hn);
      hypergraph.setNodeWeight(representative,
        hypergraph.nodeWeight(representative) - hypergraph.nodeWeight(hn));
      hypergraph.setNodePart(hn, hypergraph.partID(representative));
    }
  }

 private:
  parallel::scalable_vector<HyperedgeID> _removed_hes;
  parallel::scalable_vector<HypernodeID> _removed_hns;
  parallel::scalable_vector<HypernodeID> _representative_hns;
};
}  // namespace mt_kahypar

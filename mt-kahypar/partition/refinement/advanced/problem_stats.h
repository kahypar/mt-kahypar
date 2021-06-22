/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/range.h"

namespace mt_kahypar {

class ProblemStats {

  static constexpr size_t INVALID_IDX = std::numeric_limits<size_t>::max();

  using BlockIterator = typename vec<PartitionID>::const_iterator;

 public:
  explicit ProblemStats(const HyperedgeID num_edges,
                                const PartitionID k) :
    _k(k),
    _num_nodes_in_block(),
    _node_weight_of_block(),
    _contained_blocks(),
    _block_to_idx(k, INVALID_IDX),
    _locked_blocks(k, false),
    _num_nodes(0),
    _total_weight(0),
    _num_edges(0),
    _num_pins(0),
    _visited_hes(num_edges, false) { }

  ProblemStats(const ProblemStats&) = delete;
  ProblemStats(ProblemStats&&) = delete;

  ProblemStats & operator= (const ProblemStats &) = delete;
  ProblemStats & operator= (ProblemStats &&) = delete;

  // ####################### Accessors ######################

  HypernodeID numNodes() const {
    return _num_nodes;
  }

  HypernodeWeight totalWeight() const {
    return _total_weight;
  }

  HyperedgeID numEdges() const {
    return _num_edges;
  }

  HypernodeID numPins() const {
    return _num_pins;
  }

  PartitionID block(const PartitionID idx) const {
    ASSERT(static_cast<size_t>(idx) < _contained_blocks.size());
    return _contained_blocks[idx];
  }

  bool isBlockContained(const PartitionID block) const {
    return _block_to_idx[block] != INVALID_IDX;
  }

  PartitionID numContainedBlocks() const {
    return _contained_blocks.size();
  }

  IteratorRange<BlockIterator> containedBlocks() const {
    return IteratorRange<BlockIterator>(
      _contained_blocks.cbegin(), _contained_blocks.cend());
  }

  bool isLocked(const PartitionID block) const {
    ASSERT(block > kInvalidPartition && block < _k);
    return _locked_blocks[block];
  }

  HypernodeID numberOfNodesInBlock(const PartitionID block) const {
    ASSERT(block > kInvalidPartition && block < _k);
    const size_t idx = _block_to_idx[block];
    return idx == INVALID_IDX ? 0 : _num_nodes_in_block[idx];
  }

  HypernodeWeight nodeWeightOfBlock(const PartitionID block) const {
    ASSERT(block > kInvalidPartition && block < _k);
    const size_t idx = _block_to_idx[block];
    return idx == INVALID_IDX ? 0 : _node_weight_of_block[idx];
  }

  // ####################### Modifiers ######################

  void addBlock(const PartitionID block) {
    size_t idx = _block_to_idx[block];
    if ( idx == INVALID_IDX ) {
      idx = _contained_blocks.size();
      _block_to_idx[block] = idx;
      _contained_blocks.push_back(block);
      _num_nodes_in_block.push_back(0);
      _node_weight_of_block.push_back(0);
    }
  }

  void addNode(const HypernodeID hn,
               const PartitionID block,
               const PartitionedHypergraph& phg) {
    ASSERT(!_locked_blocks[block]);
    // Increment stats
    size_t idx = _block_to_idx[block];
    ASSERT(idx != INVALID_IDX && _contained_blocks[idx] == block);
    const HypernodeWeight node_weight = phg.nodeWeight(hn);
    ++_num_nodes_in_block[idx];
    _node_weight_of_block[idx] += node_weight;
    ++_num_nodes;
    _total_weight += node_weight;
    _num_pins += phg.nodeDegree(hn);
  }

  void addEdge(const HyperedgeID he) {
    if ( !_visited_hes[he] ) {
      ++_num_edges;
      _visited_hes[he] = true;
    }
  }

  void lockBlock(const PartitionID block) {
    ASSERT(block > kInvalidPartition && block < _k);
    _locked_blocks[block] = true;
  }

  void reset() {
    _num_nodes_in_block.clear();
    _node_weight_of_block.clear();
    _contained_blocks.clear();
    _block_to_idx.assign(_k, INVALID_IDX);
    _locked_blocks.assign(_k, false);
    _num_nodes = 0;
    _total_weight = 0;
    _num_edges = 0;
    _num_pins = 0;
    std::fill(_visited_hes.begin(), _visited_hes.end(), false);
  }

 private:
  const PartitionID _k;
  vec<HypernodeID> _num_nodes_in_block;
  vec<HypernodeWeight> _node_weight_of_block;
  vec<PartitionID> _contained_blocks;
  vec<size_t> _block_to_idx;
  vec<bool> _locked_blocks;
  HypernodeID _num_nodes;
  HypernodeWeight _total_weight;
  HyperedgeID _num_edges;
  HypernodeID _num_pins;

 public:
  vec<bool> _visited_hes;
};

inline std::ostream& operator<<(std::ostream& out, const ProblemStats& stats) {
  out << "[Nodes=" << stats.numNodes()
      << ", Edges=" << stats.numEdges()
      << ", Pins=" << stats.numPins()
      << ", Blocks=( ";
  for ( const PartitionID block : stats.containedBlocks() ) {
    out << block << " ";
  }
  out << "), Weights=( ";
  for ( const PartitionID block : stats.containedBlocks() ) {
    out << stats.nodeWeightOfBlock(block) << " ";
  }
  out << ")]";
  return out;
}

}  // namespace kahypar

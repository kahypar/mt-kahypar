/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2022 Nikolai Maas <nikolai.maas@student.kit.edu>
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

#pragma once

#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/separated_nodes.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace ds {

class SepNodesTracker {
 public:
  class Bucket {
   public:
    struct Entry {
      HypernodeID id;
      HypernodeWeight weight;
      HyperedgeWeight incident_weight;

      double density() const {
        return static_cast<double>(incident_weight) / static_cast<double>(weight);
      }

      bool operator<(const Entry& other) const {
        return density() < other.density();
      }

      bool equalStats(const Entry& other) const {
        return weight == other.weight && incident_weight == other.incident_weight;
      }
    };

    using Handle = size_t;
    static constexpr Handle empty_handle = std::numeric_limits<Handle>::max();

    explicit Bucket();

    void addNode(HypernodeID id, HypernodeWeight node_weight, HyperedgeWeight incident_weight, Handle& handle);

    void removeNode(Handle handle, Array<std::pair<PartitionID, Handle>>& handles);

    void updateWeight(HypernodeWeight new_weight);

    void reorder(Array<std::pair<PartitionID, Handle>>& handles);

    void calculateDeltasForNodeRemoval(const vec<HypernodeWeight>& node_weights,
                                       const vec<Entry>& new_nodes,
                                       vec<double>& out,
                                       Array<std::pair<PartitionID, Handle>>& handles);

    void calculateDeltasForAddingNodes(const vec<HypernodeWeight>& node_weights,
                                       const vec<Entry>& removed_nodes,
                                       vec<double>& out,
                                       Array<std::pair<PartitionID, Handle>>& handles);

    void clear();

    // ! only for testing
    bool verifyHandle(HypernodeID id, Handle handle) const  {
      return handle < _nodes.size() && _nodes[handle].id == id;
    }

    // ! only for testing
    HypernodeWeight currentWeight() const {
      return _current_weight;
    }

    // ! only for testing
    HypernodeWeight totalWeight() const {
      return _total_weight;
    }

    // ! only for testing
    size_t firstRemoved() const {
      return _first_removed;
    }

   private:
    void ensureOrdered(Array<std::pair<PartitionID, Handle>>& handles);

   private:
    vec<Entry> _nodes;
    HypernodeWeight _current_weight;
    HypernodeWeight _total_weight;
    HypernodeWeight _allowed_weight;
    size_t _first_removed;
    bool _ordered;
  };

  using BucketT = Bucket;
  using Handle = typename BucketT::Handle;

 public:
  explicit SepNodesTracker(const SeparatedNodes& s_nodes, const std::vector<HypernodeWeight>& max_part_weights, PartitionID k);

  // return whether moves need be updated globally
  bool applyMove(const SeparatedNodes& s_nodes, const vec<CAtomic<HypernodeWeight>>& part_weights,
                 const Array<CAtomic<PartitionID>>& part_ids, HypernodeID node);

  void reset();

 private:
  PartitionID currentPart(HypernodeID s_node) const {
    return _handles[s_node].first;
  }

  void setPart(HypernodeID s_node, PartitionID part) {
    _handles[s_node].first = part;
  }

  Handle& handle(HypernodeID s_node) {
    return _handles[s_node].second;
  }

  const Handle& handle(HypernodeID s_node) const {
    return _handles[s_node].second;
  }

 private:
  Array<BucketT> _buckets;
  const std::vector<HypernodeWeight> _max_part_weights;
  Array<std::pair<PartitionID, Handle>> _handles;
  HypernodeWeight _max_node_weight;
};

} // namespace ds
} // namespace mt_kahypar

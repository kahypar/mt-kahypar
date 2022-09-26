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

#include "sep_nodes_tracker.h"

namespace mt_kahypar {
namespace ds {
using Bucket = SepNodesTracker::Bucket;

Bucket::Bucket() : _nodes(), _current_weight(0), _total_weight(0),
                   _allowed_weight(0), _first_removed(0), _ordered(true) { }

void Bucket::addNode(HypernodeID id, HypernodeWeight node_weight, HyperedgeWeight incident_weight, Handle& handle) {
  handle = _nodes.size();
  _nodes.push_back({id, node_weight, incident_weight});
  _total_weight += node_weight;
  _ordered = false;
}

void Bucket::removeNode(Handle handle, Array<std::pair<PartitionID, Handle>>& handles) {
  ASSERT(handle < _nodes.size());
  _total_weight -= _nodes[handle].weight;
  std::swap(_nodes[handle], _nodes.back());
  handles[_nodes[handle].id].second = handle; // must be before pop_back()!
  _nodes.pop_back();
  _ordered = false;
}

void Bucket::updateWeight(HypernodeWeight new_weight) {
  if (new_weight != _allowed_weight) {
    _ordered = false;
    _allowed_weight = std::max(new_weight, 0);
  }
}

void Bucket::reorder(Array<std::pair<PartitionID, Handle>>& handles) {
  std::sort(_nodes.rbegin(), _nodes.rend()); // sort descending
  for (size_t i = 0; i < _nodes.size(); ++i) {
    handles[_nodes[i].id].second = i;
  }

  _current_weight = _total_weight;
  _first_removed = _nodes.size();
  while (_current_weight > _allowed_weight) {
    --_first_removed;
    _current_weight -= _nodes[_first_removed].weight;
    ASSERT(_first_removed > 0 || _current_weight == 0);
  }
  _ordered = true;
}

void Bucket::ensureOrdered(Array<std::pair<PartitionID, Handle>>& handles) {
  if (!_ordered) {
    reorder(handles);
  }
}

void Bucket::calculateDeltasForNodeRemoval(const vec<Entry>& node_weights,
                                           const vec<Entry>& new_nodes,
                                           vec<double>& out,
                                           Array<std::pair<PartitionID, Handle>>& handles) {
  ASSERT([&]{
    HypernodeWeight r_sum = 0;
    for (size_t i = 0; i < node_weights.size(); ++i) {
      if (i > 0 && node_weights[i - 1] < node_weights[i]) {
        return false;
      }
      r_sum += node_weights[i].weight;
    }
    HypernodeWeight n_sum = 0;
    for (size_t i = 0; i < new_nodes.size(); ++i) {
      if (i > 0 && new_nodes[i - 1] < new_nodes[i]) {
        return false;
      }
      n_sum += new_nodes[i].weight;
    }
    return n_sum <= r_sum;
  }());

  if (_total_weight > _allowed_weight) {
    ensureOrdered(handles);
  } else {
    _first_removed = _nodes.size();
    _current_weight = _total_weight;
  }
  out.clear();

  size_t i_left = _first_removed;
  size_t i_right = 0;
  auto next_node = [&] {
    if (i_left == _nodes.size()) {
      if (i_right == new_nodes.size()) {
        return Entry{0, std::numeric_limits<HypernodeWeight>::max(), 0}; // sentinel
      }
      return new_nodes[i_right++];
    } else if (i_right == new_nodes.size() || new_nodes[i_right] < _nodes[i_left]) {
      return _nodes[i_left++];
    }
    return new_nodes[i_right++];
  };

  Entry current_node = next_node();
  double current_density = current_node.density();
  current_node.weight -= std::max(_allowed_weight - _current_weight, 0);
  for (const Entry& e: node_weights) {
    double delta = 0;
    HypernodeWeight remaining = e.weight;
    while (remaining > 0) {
      const HypernodeWeight diff = std::min(current_node.weight, remaining);
      current_node.weight -= diff;
      remaining -= diff;
      delta += diff * current_density;

      if (current_node.weight == 0) {
        current_node = next_node();
        current_density = current_node.density();
      }
    }
    out.push_back(delta);
  }
}

void Bucket::calculateDeltasForAddingNodes(const vec<Entry>& node_weights,
                                           const vec<Entry>& removed_nodes,
                                           vec<double>& out,
                                           Array<std::pair<PartitionID, Handle>>& handles) {
  ASSERT([&]{
    for (size_t i = 0; i < node_weights.size(); ++i) {
      if (i > 0 && node_weights[i - 1] < node_weights[i]) {
        return false;
      }
    }
    for (size_t i = 0; i < removed_nodes.size(); ++i) {
      if (i > 0 && removed_nodes[i] < removed_nodes[i - 1]) {
        return false;
      }
    }
    return true;
  }());

  HypernodeWeight total_added_weight = 0;
  HypernodeWeight total_removed_weight = 0;
  for (const Entry& e: node_weights) {
    total_added_weight += e.weight;
  }
  ASSERT(total_added_weight <= _allowed_weight);
  for (const Entry& e: removed_nodes) {
    total_removed_weight += e.weight;
  }

  if (_total_weight + total_added_weight - total_removed_weight > _allowed_weight) {
    ensureOrdered(handles);
  } else {
    _first_removed = _nodes.size();
    _current_weight = _total_weight;
  }
  out.clear();

  // find starting position
  HypernodeWeight remaining_weight_upward = total_removed_weight;
  size_t index = _first_removed;
  if (index < _nodes.size()) {
    HypernodeWeight current_w = _nodes[index].weight;
    current_w -= (_allowed_weight - _current_weight);
    while (current_w < remaining_weight_upward) {
      if (current_w < remaining_weight_upward) {
        remaining_weight_upward -= current_w;
        ++index;
        current_w = (index == _nodes.size()) ? remaining_weight_upward : _nodes[index].weight;
      }
    }
  } else {
    remaining_weight_upward += (_allowed_weight - _current_weight);
  }

  size_t i_right = 0;
  auto next_node = [&] {
    ASSERT(index > 0);
    --index;
    while (i_right < removed_nodes.size()
           && _nodes[index].equalStats(removed_nodes[i_right])) {
      ASSERT(index > 0);
      --index;
      ++i_right;
    }
    ASSERT(i_right == removed_nodes.size() || _nodes[index] < removed_nodes[i_right]);
    return _nodes[index];
  };

  Entry current_node = (index == _nodes.size()) ? Entry{0, 1, 0} : _nodes[index];
  double current_density = current_node.density();
  current_node.weight = remaining_weight_upward;
  for (const Entry& e: node_weights) {
    double delta = 0;
    HypernodeWeight remaining = e.weight;
    while (remaining > 0) {
      if (current_node.weight == 0) {
        current_node = next_node();
        current_density = current_node.density();
      }

      const HypernodeWeight diff = std::min(current_node.weight, remaining);
      current_node.weight -= diff;
      remaining -= diff;
      delta += diff * current_density;
    }
    out.push_back(delta);
  }
}

void Bucket::clear() {
  _nodes.clear();
  _current_weight = 0;
  _total_weight = 0;
  _allowed_weight = 0;
  _first_removed = 0;
}


SepNodesTracker::SepNodesTracker(const SeparatedNodes& s_nodes,
                                 const std::vector<HypernodeWeight>& max_part_weights,
                                 PartitionID k) :
        _buckets(),
        _handles() {
        _max_part_weights(&max_part_weights),
  ASSERT(max_part_weights.size() == static_cast<size_t>(k));
  _buckets.resize(k, Bucket(), false);
  _handles.resize(s_nodes.numNodes(), {kInvalidPartition, BucketT::empty_handle}, false);
}

bool SepNodesTracker::applyMove(const SeparatedNodes& s_nodes, const vec<CAtomic<HypernodeWeight>>& part_weights,
                                const Array<CAtomic<PartitionID>>& part_ids, HypernodeID node) {
  ASSERT(part_weights.size() == _buckets.size());
  // update bucket state regarding separated nodes
  for (const SeparatedNodes::Edge& outward: s_nodes.outwardEdges(node)) {
    const HypernodeID u = outward.target;
    HyperedgeWeight total_incident_weight = 0;
    HyperedgeWeight invalid_incident_weight = 0;
    vec<HyperedgeWeight> block_scores(part_weights.size(), 0);
    for (const SeparatedNodes::Edge& e: s_nodes.inwardEdges(u)) {
      total_incident_weight += e.weight;
      PartitionID part = part_ids[e.target];
      if (part == kInvalidPartition) {
        invalid_incident_weight += e.weight;
      } else {
        ASSERT(static_cast<size_t>(part) < block_scores.size());
        block_scores[part] += e.weight;
      }
    }

    PartitionID max_part = kInvalidPartition;
    HyperedgeWeight max_score = 0;
    HyperedgeWeight second_score = 0;
    for (size_t part = 0; part < block_scores.size(); ++part) {
      if (block_scores[part] >= max_score) {
        second_score = max_score;
        max_score = block_scores[part];
        max_part = part;
      } else if (block_scores[part] > second_score) {
        second_score = block_scores[part];
      }
    }
    ASSERT(max_score >= second_score);

    const HyperedgeWeight value = max_score - second_score;
    const PartitionID target_part = (value >= invalid_incident_weight) ? max_part : kInvalidPartition;

    if (currentPart(u) != kInvalidPartition) {
      _buckets[currentPart(u)].removeNode(handle(u), _handles);
    }
    if (target_part != kInvalidPartition) {
      _buckets[target_part].addNode(u, s_nodes.nodeWeight(u), value, handle(u));
    }
    setPart(u, target_part);
  }

  for (size_t part = 0; part < _buckets.size(); ++part) {
    _buckets[part].updateWeight((*_max_part_weights)[part] - part_weights[part].load());
  }
  return true; // TODO
}

void SepNodesTracker::reset() {
  for (Bucket& bucket: _buckets) {
    bucket.clear();
  }
  for (auto& entry: _handles) {
    entry = {kInvalidPartition, BucketT::empty_handle};
  }
}

} // namespace ds
} // namespace mt_kahypar

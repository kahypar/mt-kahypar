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

#include <cmath>

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
      Entry result = _nodes[i_left];
      // if (i_left == _first_removed) {
      //   // TODO
      //   result.weight -= std::max(_allowed_weight - _current_weight, 0);
      // }
      ++i_left;
      return result;
    }
    return new_nodes[i_right++];
  };

  Entry current_node = next_node();
  double current_density = current_node.density();
  // current_node.weight -= std::max(_allowed_weight - _current_weight, 0);
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
    if (index == 0) {
      return Entry{0, 10, std::numeric_limits<HyperedgeWeight>::max()}; // sentinel
    }
    --index;
    while (i_right < removed_nodes.size()
           && _nodes[index].equalStats(removed_nodes[i_right])) {
      if (index == 0) {
        return Entry{0, 10, std::numeric_limits<HyperedgeWeight>::max()}; // sentinel
      }
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
        _max_part_weights(&max_part_weights),
        _handles(),
        _removed_from(), _inserted_from(),
        _removed_to(), _inserted_to() {
  ASSERT(max_part_weights.size() == static_cast<size_t>(k));
  _buckets.resize(k, Bucket(), false);
  _handles.resize(s_nodes.numNodes(), {kInvalidPartition, BucketT::empty_handle}, false);
}

bool SepNodesTracker::applyMove(const SeparatedNodes& s_nodes, const vec<CAtomic<HypernodeWeight>>& part_weights,
                                const Array<CAtomic<PartitionID>>& part_ids, HypernodeID node) {
  ASSERT(part_weights.size() == _buckets.size());
  // update bucket state regarding separated nodes
  Array<HyperedgeWeight> block_scores;
  for (const SeparatedNodes::Edge& outward: s_nodes.outwardEdges(node)) {
    const HypernodeID u = outward.target;
    auto [target_part, value] = calculateNewPartOfSeparated(u, s_nodes, part_ids, part_weights.size(), block_scores);
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


HyperedgeWeight SepNodesTracker::rateMove(const SeparatedNodes& s_nodes, const Array<CAtomic<PartitionID>>& part_ids, PartitionID k,
                                          HypernodeID node, HypernodeWeight weight, PartitionID from, PartitionID to) {
  // TODO: k > 2 (removing from other blocks)
  ASSERT(k == 2);
  _removed_from.clear();
  _inserted_from.clear();
  _removed_to.clear();
  _inserted_to.clear();

  HyperedgeWeight moved_benefit = 0;
  // collect added and removed nodes
  Entry dummy{kInvalidHypernode, weight, std::numeric_limits<HyperedgeWeight>::max()};
  _removed_from.push_back(dummy);
  _inserted_to.push_back(dummy);
  Array<HyperedgeWeight> block_scores;
  for (const SeparatedNodes::Edge& outward: s_nodes.outwardEdges(node)) {
    const HypernodeID u = outward.target;
    if (currentPart(u) != kInvalidPartition) {
      Entry old_entry{u, s_nodes.nodeWeight(u), _buckets[currentPart(u)].get(handle(u)).incident_weight};
      if (currentPart(u) == from) {
        _removed_from.push_back(old_entry);
      } else {
        ASSERT(currentPart(u) == to);
        _removed_to.push_back(old_entry);
        if (_buckets[to].isIncluded(handle(u), _handles)) {
          moved_benefit -= old_entry.incident_weight;
        } else {
          moved_benefit += old_entry.incident_weight;
        }
      }
    }

    auto [target_part, value] = calculateNewPartOfSeparated(u, s_nodes, part_ids, k, block_scores,
                                                            from, to, outward.weight);
    if (target_part != kInvalidPartition) {
      Entry new_entry{u, s_nodes.nodeWeight(u), value};
      if (target_part == from) {
        _inserted_from.push_back(new_entry);
        moved_benefit -= new_entry.incident_weight;
      } else {
        ASSERT(target_part == to);
        _inserted_to.push_back(new_entry);
      }
    }
  }

  std::sort(_removed_from.rbegin(), _removed_from.rend());
  std::sort(_inserted_from.rbegin(), _inserted_from.rend());
  std::sort(_inserted_to.rbegin(), _inserted_to.rend());
  std::sort(_removed_to.begin(), _removed_to.end()); // note: needs ascending order

  double result = 0;
  // calculate from part
  _buckets[from].calculateDeltasForNodeRemoval(_removed_from, _inserted_from, _out, _handles);
  ASSERT(_removed_from.size() == _out.size());
  bool not_in_bucket = false;
  for (size_t i = 0; i < _out.size(); ++i) {
    if (not_in_bucket || static_cast<double>(_removed_from[i].incident_weight) < _out[i]) {
      not_in_bucket = true;
      moved_benefit += _removed_from[i].incident_weight;
    } else {
      moved_benefit -= (i == 0) ? 0 : _removed_from[i].incident_weight;
      result += _out[i];
    }
  }
  // calculate to part
  ASSERT(weight <= _buckets[to].allowedWeight());
  _buckets[to].calculateDeltasForAddingNodes(_inserted_to, _removed_to, _out, _handles);
  ASSERT(_inserted_to.size() == _out.size());
  for (size_t i = 0; i < _out.size(); ++i) {
    if (static_cast<double>(_inserted_to[i].incident_weight) < _out[i]) {
      moved_benefit -= _inserted_to[i].incident_weight;
    } else {
      moved_benefit += (i == 0) ? 0 : _inserted_to[i].incident_weight;
      result -= _out[i];
    }
  }
  ASSERT(moved_benefit % 2 == 0);
  return static_cast<HyperedgeWeight>(std::round(result)) + moved_benefit / 2;
}

void SepNodesTracker::reset() {
  for (Bucket& bucket: _buckets) {
    bucket.clear();
  }
  for (auto& entry: _handles) {
    entry = {kInvalidPartition, BucketT::empty_handle};
  }
}

std::pair<PartitionID, HyperedgeWeight>
SepNodesTracker::calculateNewPartOfSeparated(HypernodeID separated, const SeparatedNodes& s_nodes,
                                             const Array<CAtomic<PartitionID>>& part_ids,
                                             PartitionID k, Array<HyperedgeWeight>& block_scores) {
  return calculateNewPartOfSeparated(separated, s_nodes, part_ids, k, block_scores,
                                     kInvalidPartition, kInvalidPartition, 0);
}

std::pair<PartitionID, HyperedgeWeight>
SepNodesTracker::calculateNewPartOfSeparated(HypernodeID separated, const SeparatedNodes& s_nodes,
                                             const Array<CAtomic<PartitionID>>& part_ids,
                                             PartitionID k, Array<HyperedgeWeight>& block_scores,
                                             PartitionID from, PartitionID to, HyperedgeWeight diff) {
  HyperedgeWeight total_incident_weight = 0;
  HyperedgeWeight invalid_incident_weight = 0;
  block_scores.assign(k, 0, false);
  for (const SeparatedNodes::Edge& e: s_nodes.inwardEdges(separated)) {
    total_incident_weight += e.weight;
    PartitionID part = part_ids[e.target];
    if (part == kInvalidPartition) {
      invalid_incident_weight += e.weight;
    } else {
      ASSERT(static_cast<size_t>(part) < block_scores.size());
      block_scores[part] += e.weight;
    }
  }
  if (from != kInvalidPartition) {
    block_scores[from] -= diff;
  } else {
    invalid_incident_weight -= diff;
  }
  if (to != kInvalidPartition) {
    block_scores[to] += diff;
  } else {
    invalid_incident_weight += diff;
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
  return {target_part, value};
}

} // namespace ds
} // namespace mt_kahypar

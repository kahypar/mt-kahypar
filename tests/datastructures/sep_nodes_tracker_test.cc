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

#include "gmock/gmock.h"

#include <algorithm>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/separated_nodes.h"
#include "mt-kahypar/datastructures/sep_nodes_tracker.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

using Handle = SepNodesTracker::Handle;
using Bucket = SepNodesTracker::BucketT;
using Entry = SepNodesTracker::Entry;

class ABucket: public Test {
 public:
  void setupForBucket(HypernodeID num_nodes) {
    handles = Array<std::pair<PartitionID, Handle>>();
    handles.resize(num_nodes, {kInvalidPartition, Bucket::empty_handle});
    bucket = Bucket();
  }

  void addNode(HypernodeID id, HypernodeWeight node_weight, HyperedgeWeight incident_weight) {
    bucket.addNode(id, node_weight, incident_weight, handles[id].second);
    handles[id].first = 0;
  }

  void removeNode(HypernodeID id) {
    bucket.removeNode(handles[id].second, handles);
    handles[id].first = kInvalidPartition;
  }

  void setNodesForBucket(const vec<std::pair<HypernodeWeight, HyperedgeWeight>>& nodes) {
    handles = Array<std::pair<PartitionID, Handle>>();
    handles.resize(nodes.size(), {0, Bucket::empty_handle});
    bucket = Bucket();
    for (size_t i = 0; i < nodes.size(); ++i) {
      bucket.addNode(i, nodes[i].first, nodes[i].second, handles[i].second);
    }
  }

  void verifyHandles() {
    for (size_t i = 0; i < handles.size(); ++i) {
      if (handles[i].first == 0) {
        ASSERT_TRUE(bucket.verifyHandle(i, handles[i].second));
      }
    }
  }

  Bucket bucket;
  Array<std::pair<PartitionID, Handle>> handles;
};

class ATracker: public Test {
 public:
  void setupWithUnitWeights(HypernodeID num_graph_nodes,
                            vec<vec<std::pair<HypernodeID, HyperedgeWeight>>> separated,
                            const std::vector<HypernodeWeight>& max_part_weights,
                            const vec<PartitionID>& input_part_ids,
                            const Array<HypernodeWeight>& node_weights = {}) {
    ASSERT(num_graph_nodes == input_part_ids.size());
    Hypergraph graph = HypergraphFactory::construct(num_graph_nodes, 0, {});
    vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> nodes;
    vec<SeparatedNodes::Edge> s_edges;
    for (const auto& edges_of_node: separated) {
      nodes.emplace_back(kInvalidHypernode, s_edges.size(), 1);
      for (const auto& [target, weight]: edges_of_node) {
        s_edges.emplace_back(target, weight);
      }
    }
    s_nodes = SeparatedNodes(num_graph_nodes);
    s_nodes.addNodes(graph, nodes, s_edges);
    s_nodes.revealAll();
    s_nodes.initializeOutwardEdges();

    tracker = SepNodesTracker(s_nodes, max_part_weights, max_part_weights.size());
    vec<CAtomic<HypernodeWeight>> part_weights;
    part_weights.resize(max_part_weights.size(), CAtomic<HypernodeWeight>(0));
    part_ids = Array<CAtomic<PartitionID>>();
    part_ids.resize(input_part_ids.size(), CAtomic<PartitionID>(kInvalidPartition));
    for (HypernodeID node = 0; node < num_graph_nodes; ++node) {
      part_ids[node].store(input_part_ids[node]);
      if (!node_weights.empty()) {
        part_weights[input_part_ids[node]] += node_weights[node];
      }
      tracker.applyMove(s_nodes, part_weights, part_ids, node);
    }
  }

 public:
  SeparatedNodes s_nodes;
  SepNodesTracker tracker;
  Array<CAtomic<PartitionID>> part_ids;
};

vec<Entry> entries(vec<HypernodeWeight>&& weights) {
  vec<Entry> result;
  for (const HypernodeWeight& w: weights) {
    result.push_back(Entry{kInvalidHypernode, w, w});
  }
  return result;
}

TEST_F(ABucket, AddsNodesAndReorders1) {
  setupForBucket(3);
  addNode(0, 1, 3);
  addNode(1, 1, 2);
  addNode(2, 1, 1);
  ASSERT_EQ(3, bucket.totalWeight());
  verifyHandles();

  bucket.updateWeight(3);
  bucket.reorder(handles);
  ASSERT_EQ(3, bucket.currentWeight());
  ASSERT_EQ(3, bucket.totalWeight());
  ASSERT_EQ(3, bucket.firstRemoved());
  verifyHandles();
}

TEST_F(ABucket, AddsNodesAndReorders2) {
  setupForBucket(4);
  addNode(0, 1, 3);
  addNode(1, 2, 2);
  addNode(2, 2, 1);
  addNode(3, 1, 2);
  ASSERT_EQ(6, bucket.totalWeight());
  verifyHandles();

  bucket.updateWeight(3);
  bucket.reorder(handles);
  ASSERT_EQ(2, bucket.currentWeight());
  ASSERT_EQ(6, bucket.totalWeight());
  ASSERT_EQ(2, bucket.firstRemoved());
  verifyHandles();
}

TEST_F(ABucket, AddRemoveReorder) {
  setupForBucket(4);
  addNode(0, 2, 3);
  addNode(1, 4, 2);
  addNode(2, 1, 1);
  addNode(3, 5, 3);
  ASSERT_EQ(12, bucket.totalWeight());
  removeNode(1);
  removeNode(2);
  verifyHandles();
  ASSERT_EQ(7, bucket.totalWeight());
  bucket.updateWeight(2);
  bucket.reorder(handles);
  ASSERT_EQ(2, bucket.currentWeight());
  ASSERT_EQ(1, bucket.firstRemoved());

  addNode(1, 2, 7);
  addNode(2, 1, 2);
  removeNode(0);
  removeNode(3);
  verifyHandles();
  bucket.updateWeight(1);
  bucket.reorder(handles);
  verifyHandles();
  ASSERT_EQ(0, bucket.currentWeight());
  ASSERT_EQ(0, bucket.firstRemoved());
}

TEST_F(ABucket, DeltaForRemoval) {
  setNodesForBucket({ {1, 1}, {1, 2}, {2, 5}, {1, 3} });
  bucket.updateWeight(5);
  vec<double> result;

  // empty
  bucket.calculateDeltasForNodeRemoval({}, {}, result, handles);
  ASSERT_EQ(0, result.size());
  bucket.calculateDeltasForNodeRemoval(entries({1}), {}, result, handles);
  ASSERT_EQ(1, result.size());
  ASSERT_EQ(0.0, result[0]);
  bucket.calculateDeltasForNodeRemoval(entries({1, 2, 3}), {}, result, handles);
  ASSERT_EQ(3, result.size());
  ASSERT_EQ(0.0, result[0]);
  ASSERT_EQ(0.0, result[1]);
  ASSERT_EQ(0.0, result[2]);

  // some new nodes
  bucket.calculateDeltasForNodeRemoval(entries({3}), {Entry{0, 1, 2}, Entry{0, 1, 1}}, result, handles);
  ASSERT_EQ(1, result.size());
  ASSERT_EQ(3.0, result[0]);
  bucket.calculateDeltasForNodeRemoval(entries({1, 1}), {Entry{0, 1, 2}, Entry{0, 1, 1}}, result, handles);
  ASSERT_EQ(2, result.size());
  ASSERT_EQ(2.0, result[0]);
  ASSERT_EQ(1.0, result[1]);

  // now use a full one
  bucket.updateWeight(3);
  bucket.calculateDeltasForNodeRemoval(entries({1}), {}, result, handles);
  ASSERT_EQ(1, result.size());
  ASSERT_EQ(2.0, result[0]);
  bucket.calculateDeltasForNodeRemoval(entries({1, 2, 1}), {}, result, handles);
  ASSERT_EQ(3, result.size());
  ASSERT_EQ(2.0, result[0]);
  ASSERT_EQ(1.0, result[1]);
  ASSERT_EQ(0.0, result[2]);

  // combined
  bucket.calculateDeltasForNodeRemoval(entries({2, 1, 1}), {Entry{0, 2, 3}}, result, handles);
  ASSERT_EQ(3, result.size());
  ASSERT_EQ(3.5, result[0]);
  ASSERT_EQ(1.5, result[1]);
  ASSERT_EQ(1.0, result[2]);
  bucket.updateWeight(7);
  bucket.calculateDeltasForNodeRemoval(entries({2}), {Entry{0, 2, 2}}, result, handles);
  ASSERT_EQ(1, result.size());
  ASSERT_EQ(2, result[0]);
  bucket.updateWeight(2);
  bucket.calculateDeltasForNodeRemoval(entries({2}), {}, result, handles);
  ASSERT_EQ(1, result.size());
  ASSERT_EQ(4.5, result[0]);
}

TEST_F(ABucket, DeltaForAdding) {
  setNodesForBucket({ {1, 1}, {2, 5}, {1, 3} });
  bucket.updateWeight(7);
  vec<double> result;

  // empty
  bucket.calculateDeltasForAddingNodes({}, {}, result, handles);
  ASSERT_EQ(0, result.size());
  bucket.calculateDeltasForAddingNodes(entries({2, 1}), {}, result, handles);
  ASSERT_EQ(2, result.size());
  ASSERT_EQ(0.0, result[0]);
  ASSERT_EQ(0.0, result[1]);
  bucket.calculateDeltasForAddingNodes(entries({2, 3}), {Entry{0, 2, 5}}, result, handles);
  ASSERT_EQ(2, result.size());
  ASSERT_EQ(0.0, result[0]);
  ASSERT_EQ(0.0, result[1]);

  // now it gets interesting
  bucket.updateWeight(5);
  bucket.calculateDeltasForAddingNodes(entries({2, 1}), {}, result, handles);
  ASSERT_EQ(2, result.size());
  ASSERT_EQ(1.0, result[0]);
  ASSERT_EQ(2.5, result[1]);
  bucket.calculateDeltasForAddingNodes(entries({1, 2, 2}), {}, result, handles);
  ASSERT_EQ(3, result.size());
  ASSERT_EQ(0.0, result[0]);
  ASSERT_EQ(3.5, result[1]);
  ASSERT_EQ(5.5, result[2]);
  bucket.calculateDeltasForAddingNodes(entries({3, 2}), {Entry{0, 1, 1}, Entry{0, 1, 3}}, result, handles);
  ASSERT_EQ(2, result.size());
  ASSERT_EQ(0.0, result[0]);
  ASSERT_EQ(5.0, result[1]);

  setNodesForBucket({ {1, 1}, {2, 4}, {2, 5}, {1, 3} });
  bucket.updateWeight(4);
  bucket.calculateDeltasForAddingNodes(entries({1, 2}), {Entry{0, 1, 1}, Entry{0, 2, 4}, Entry{0, 1, 3}}, result, handles);
  ASSERT_EQ(2, result.size());
  ASSERT_EQ(0.0, result[0]);
  ASSERT_EQ(2.5, result[1]);
}

TEST_F(ATracker, InitializesWithCorrectStats) {
  setupWithUnitWeights(1, { {{0, 1}} }, {2, 2}, {0});
  ASSERT_EQ(1, tracker.buckets()[0].totalWeight());
  ASSERT_EQ(0, tracker.buckets()[1].totalWeight());
  setupWithUnitWeights(2, { {{0, 1}, {1, 2}}, {{0, 2}, {1, 1}} }, {2, 2}, {0, 1});
  ASSERT_EQ(1, tracker.buckets()[0].totalWeight());
  ASSERT_EQ(1, tracker.buckets()[1].totalWeight());
}

TEST_F(ATracker, UnitWeightBasicCases) {
  setupWithUnitWeights(1, { {{0, 1}} }, {2, 2}, {0});
  // ASSERT_EQ(0, tracker.rateMove(s_nodes, part_ids, 2, 0, 0, 0, 0));
  ASSERT_EQ(0, tracker.rateMove(s_nodes, part_ids, 2, 0, 0, 0, 1));
  setupWithUnitWeights(2, { {{0, 1}, {1, 1}} }, {2, 2}, {0, 1});
  ASSERT_EQ(1, tracker.rateMove(s_nodes, part_ids, 2, 0, 0, 0, 1));
  setupWithUnitWeights(2, { {{0, 1}, {1, 1}} }, {2, 2}, {0, 0});
  ASSERT_EQ(-1, tracker.rateMove(s_nodes, part_ids, 2, 0, 0, 0, 1));
  setupWithUnitWeights(2, { {{0, 1}, {1, 3}} }, {2, 2}, {0, 0});
  ASSERT_EQ(-1, tracker.rateMove(s_nodes, part_ids, 2, 0, 0, 0, 1));
  setupWithUnitWeights(2, { {{0, 1}}, {{1, 1}} }, {1, 1}, {0, 0});
  ASSERT_EQ(1, tracker.rateMove(s_nodes, part_ids, 2, 0, 0, 0, 1));
  setupWithUnitWeights(2, { {{0, 2}}, {{1, 1}}, {{1, 1}} }, {1, 1}, {0, 0});
  ASSERT_EQ(1, tracker.rateMove(s_nodes, part_ids, 2, 0, 0, 0, 1));
  setupWithUnitWeights(2, { {{0, 2}}, {{1, 1}} }, {1, 1}, {0, 1});
  ASSERT_EQ(-1, tracker.rateMove(s_nodes, part_ids, 2, 0, 0, 0, 1));
  setupWithUnitWeights(2, { {{0, 2}}, {{1, 1}}, {{1, 1}} }, {1, 1}, {0, 1});
  ASSERT_EQ(-1, tracker.rateMove(s_nodes, part_ids, 2, 0, 0, 0, 1));
  setupWithUnitWeights(2, { {{0, 1}}, {{0, 1}}, {{1, 2}} }, {2, 2}, {0, 0});
  ASSERT_EQ(1, tracker.rateMove(s_nodes, part_ids, 2, 0, 0, 0, 1));

  // with node weight
  setupWithUnitWeights(1, { {{0, 1}} }, {2, 2}, {0});
  ASSERT_EQ(0, tracker.rateMove(s_nodes, part_ids, 2, 0, 1, 0, 1));
  setupWithUnitWeights(2, { {{0, 2}}, {{1, 1}}, {{1, 1}}, {{1, 1}} }, {1, 3}, {0, 0});
  ASSERT_EQ(3, tracker.rateMove(s_nodes, part_ids, 2, 0, 2, 0, 1));
  setupWithUnitWeights(2, { {{0, 2}}, {{1, 1}}, {{1, 1}} }, {1, 2}, {0, 1});
  ASSERT_EQ(-2, tracker.rateMove(s_nodes, part_ids, 2, 0, 1, 0, 1));
}

}
} // namespace mt_kahypar

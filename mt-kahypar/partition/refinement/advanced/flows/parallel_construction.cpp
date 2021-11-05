/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/partition/refinement/advanced/flows/parallel_construction.h"

#include "tbb/concurrent_queue.h"

#include "mt-kahypar/partition/refinement/advanced/flows/flow_common.h"
#include "mt-kahypar/parallel/stl/scalable_queue.h"

namespace mt_kahypar {

ParallelConstruction::TmpHyperedge ParallelConstruction::DynamicIdenticalNetDetection::get(const size_t he_hash,
                                                                                           const vec<whfc::Node>& pins) {
  const size_t* bucket_idx = _he_hashes.get_if_contained(he_hash);
  if ( bucket_idx ) {
    // There exists already some hyperedges with the same hash
    for ( const TmpHyperedge& tmp_e : _hash_buckets[*bucket_idx] ) {
      // Check if there is some hyperedge equal to he
      if ( _flow_hg.pinCount(tmp_e.bucket, tmp_e.e) == pins.size() ) {
        bool is_identical = true;
        size_t idx = 0;
        for ( const whfc::FlowHypergraph::Pin& u : _flow_hg.tmpPinsOf(tmp_e.bucket, tmp_e.e) ) {
          if ( u.pin != pins[idx++] ) {
            is_identical = false;
            break;
          }
        }
        if ( is_identical ) {
          return tmp_e;
        }
      }
    }
  }
  return TmpHyperedge { std::numeric_limits<size_t>::max(), whfc::invalidHyperedge };
}

void ParallelConstruction::DynamicIdenticalNetDetection::add(const TmpHyperedge& tmp_he,
                                                             const size_t he_hash) {
  const size_t* bucket_idx = _he_hashes.get_if_contained(he_hash);
  // There is no hyperedge currently identical to he
  if ( bucket_idx ) {
    // If there already exist hyperedges with the same hash,
    // we insert he into the corresponding bucket.
    _hash_buckets[*bucket_idx].push_back(tmp_he);
  } else {
    // Otherwise, we create a new bucket (or reuse an existing)
    // and insert he into the hash table
    size_t idx = _used_entries++;
    ASSERT(idx < _hash_buckets.size());
    _hash_buckets[idx].clear();
    _hash_buckets[idx].push_back(tmp_he);
    _he_hashes[he_hash] = idx;
  }
}

FlowProblem ParallelConstruction::constructFlowHypergraph(const PartitionedHypergraph& phg,
                                                          const Subhypergraph& sub_hg,
                                                          const PartitionID block_0,
                                                          const PartitionID block_1,
                                                          vec<HypernodeID>& whfc_to_node) {
  ASSERT(block_0 != kInvalidPartition && block_1 != kInvalidPartition);
  FlowProblem flow_problem;
  flow_problem.total_cut = 0;
  flow_problem.non_removable_cut = 0;
  _node_to_whfc.clear();

  utils::Timer::instance().start_timer("allocate_nodes", "Allocate Nodes", true);
  tbb::parallel_invoke([&]() {
    _node_to_whfc.setMaxSize(sub_hg.numNodes());
  }, [&] {
    whfc_to_node.resize(sub_hg.numNodes() + 2);
  }, [&] {
    _flow_hg.allocateNodes(sub_hg.numNodes() + 2);
  }, [&] {
    _identical_nets.reset(sub_hg.hes.size());
  });
  utils::Timer::instance().stop_timer("allocate_nodes");

  if ( _context.refinement.advanced.flows.determine_distance_from_cut ) {
    _cut_hes.clear();
  }

  utils::Timer::instance().start_timer("add_nodes", "Add Nodes", true);
  // Add refinement nodes to flow network
  tbb::parallel_invoke([&] {
    // Add source nodes
    flow_problem.source = whfc::Node(0);
    whfc_to_node[flow_problem.source] = kInvalidHypernode;
    _flow_hg.nodeWeight(flow_problem.source) = whfc::NodeWeight(
      std::max(0, phg.partWeight(block_0) - sub_hg.weight_of_block_0));
    tbb::parallel_for(0UL, sub_hg.nodes_of_block_0.size(), [&](const size_t i) {
      const HypernodeID hn = sub_hg.nodes_of_block_0[i];
      const whfc::Node u(1 + i);
      whfc_to_node[u] = hn;
      _node_to_whfc[hn] = u;
      _flow_hg.nodeWeight(u) = whfc::NodeWeight(phg.nodeWeight(hn));
    });
  }, [&] {
    // Add sink nodes
    flow_problem.sink = whfc::Node(sub_hg.nodes_of_block_0.size() + 1);
    whfc_to_node[flow_problem.sink] = kInvalidHypernode;
    _flow_hg.nodeWeight(flow_problem.sink) = whfc::NodeWeight(
      std::max(0, phg.partWeight(block_1) - sub_hg.weight_of_block_1));
    tbb::parallel_for(0UL, sub_hg.nodes_of_block_1.size(), [&](const size_t i) {
      const HypernodeID hn = sub_hg.nodes_of_block_1[i];
      const whfc::Node u(flow_problem.sink + 1 + i);
      whfc_to_node[u] = hn;
      _node_to_whfc[hn] = u;
      _flow_hg.nodeWeight(u) = whfc::NodeWeight(phg.nodeWeight(hn));
    });
  });
  flow_problem.weight_of_block_0 = _flow_hg.nodeWeight(flow_problem.source) + sub_hg.weight_of_block_0;
  flow_problem.weight_of_block_1 = _flow_hg.nodeWeight(flow_problem.sink) + sub_hg.weight_of_block_1;
  utils::Timer::instance().stop_timer("add_nodes");

  utils::Timer::instance().start_timer("allocate_pins", "Allocate Pins", true);
  const HyperedgeID max_hyperedges = sub_hg.hes.size();
  const HypernodeID max_pins = sub_hg.num_pins + max_hyperedges;
  _flow_hg.allocateHyperedgesAndPins(max_hyperedges, max_pins);
  utils::Timer::instance().stop_timer("allocate_pins");

  // Add hyperedge to flow network and configure source and sink
  auto push_into_tmp_pins = [&](vec<whfc::Node>& tmp_pins, const whfc::Node pin,
                                size_t& current_hash, const bool is_source_or_sink) {
    tmp_pins.push_back(pin);
    current_hash += kahypar::math::hash(pin);
    if ( is_source_or_sink ) {
      // According to Lars: Adding to source or sink to the start of
      // each pin list improves running time
      std::swap(tmp_pins[0], tmp_pins.back());
    }
  };

  utils::Timer::instance().start_timer("add_pins", "Add Pins", true);
  _flow_hg.setNumCSRBuckets(NUM_CSR_BUCKETS);
  const size_t step = max_hyperedges / NUM_CSR_BUCKETS + (max_hyperedges % NUM_CSR_BUCKETS != 0);
  tbb::parallel_for(0UL, NUM_CSR_BUCKETS, [&](const size_t idx) {
    const size_t start = std::min(step * idx, static_cast<size_t>(max_hyperedges));
    const size_t end = std::min(step * (idx + 1), static_cast<size_t>(max_hyperedges));
    const size_t num_hes = end - start;
    size_t num_pins = 0;
    for ( size_t i = start; i < end; ++i ) {
      const HyperedgeID he = sub_hg.hes[i];
      num_pins += phg.edgeSize(he) + 1;
    }
    _flow_hg.initializeCSRBucket(idx, num_hes, num_pins);

    whfc::Hyperedge e(0);
    size_t pin_idx = 0;
    vec<whfc::Node>& tmp_pins = _tmp_pins.local();
    for ( size_t i = start; i < end; ++i ) {
      const HyperedgeID he = sub_hg.hes[i];
      if ( !canHyperedgeBeDropped(phg, he, block_0, block_1) ) {
        tmp_pins.clear();
        size_t he_hash = 0;
        bool connectToSource = false;
        bool connectToSink = false;
        const HyperedgeWeight he_weight = phg.edgeWeight(he);
        if ( phg.pinCountInPart(he, block_0) > 0 && phg.pinCountInPart(he, block_1) > 0 ) {
          __atomic_fetch_add(&flow_problem.total_cut, he_weight, __ATOMIC_RELAXED);
        }
        for ( const HypernodeID& pin : phg.pins(he) ) {
          whfc::Node* whfc_pin = _node_to_whfc.get_if_contained(pin);
          if ( whfc_pin ) {
            push_into_tmp_pins(tmp_pins, *whfc_pin, he_hash, false);
          } else {
            const PartitionID pin_block = phg.partID(pin);
            connectToSource |= pin_block == block_0;
            connectToSink |= pin_block == block_1;
          }
        }

        const bool empty_hyperedge = tmp_pins.size() == 0;
        const bool connected_to_source_and_sink = connectToSource && connectToSink;
        if ( connected_to_source_and_sink ) {
          // Hyperedge is connected to source and sink which means we can not remove it
          // from the cut with the current flow problem => remove he from flow problem
          __atomic_fetch_add(&flow_problem.non_removable_cut, he_weight, __ATOMIC_RELAXED);
        } else if ( !empty_hyperedge ) {
          if ( connectToSource ) {
            push_into_tmp_pins(tmp_pins, flow_problem.source, he_hash, true);
          } else if ( connectToSink ) {
            push_into_tmp_pins(tmp_pins, flow_problem.sink, he_hash, true);
          }

          // Sort pins for identical net detection
          std::sort( tmp_pins.begin() +
                   ( tmp_pins[0] == flow_problem.source ||
                     tmp_pins[0] == flow_problem.sink), tmp_pins.end());

          if ( tmp_pins.size() > 1 ) {
            const TmpHyperedge identical_net = _identical_nets.get(he_hash, tmp_pins);
            if ( identical_net.e == whfc::invalidHyperedge ) {
              const size_t pin_start = pin_idx;
              const size_t pin_end = pin_start + tmp_pins.size();
              for ( size_t i = 0; i < tmp_pins.size(); ++i ) {
                _flow_hg.addPin(tmp_pins[i], idx, pin_idx++);
              }
              TmpHyperedge tmp_e { idx, e++ };
              if ( _context.refinement.advanced.flows.determine_distance_from_cut &&
                  phg.pinCountInPart(he, block_0) > 0 && phg.pinCountInPart(he, block_1) > 0 ) {
                _cut_hes.push_back(tmp_e);
              }
              _flow_hg.finishHyperedge(tmp_e.e, he_weight, idx, pin_start, pin_end);
              _identical_nets.add(tmp_e, he_hash);
            } else {
              // Current hyperedge is identical to an already added
              __atomic_fetch_add(&_flow_hg.capacity(identical_net.bucket, identical_net.e), he_weight, __ATOMIC_RELAXED);
            }
          }
        }
      }
    }
  });

  tbb::parallel_for(0UL, NUM_CSR_BUCKETS, [&](const size_t idx) {
    _flow_hg.finalizeCSRBucket(idx);
  });
  _flow_hg.finalizeHyperedges();
  utils::Timer::instance().stop_timer("add_pins");

  if ( _flow_hg.nodeWeight(flow_problem.source) == 0 ||
       _flow_hg.nodeWeight(flow_problem.sink) == 0 ) {
    // Source or sink not connected to vertices in the flow problem
    flow_problem.non_removable_cut = 0;
    flow_problem.total_cut = 0;
  } else {
    utils::Timer::instance().start_timer("finalize", "Finalize", true);
    _flow_hg.finalizeParallel();
    utils::Timer::instance().stop_timer("finalize");

    if ( _context.refinement.advanced.flows.determine_distance_from_cut ) {
      // Determine the distance of each node contained in the flow network from the cut.
      // This technique improves piercing decision within the WHFC framework.
      utils::Timer::instance().start_timer("distance_from_cut", "Distance From Cut", true);
      determineDistanceFromCut(phg, flow_problem.source,
        flow_problem.sink, block_0, block_1, whfc_to_node);
      utils::Timer::instance().stop_timer("distance_from_cut");
    }
  }

  DBG << "Flow Hypergraph [ Nodes =" << _flow_hg.numNodes()
      << ", Edges =" << _flow_hg.numHyperedges()
      << ", Pins =" << _flow_hg.numPins()
      << ", Blocks = (" << block_0 << "," << block_1 << ") ]";

  return flow_problem;
}

namespace {
template<typename T>
class BFSQueue {

 public:
  explicit BFSQueue(const size_t num_threads) :
    _q(num_threads) { }

  bool empty() {
    bool is_empty = true;
    for ( size_t i = 0; i < _q.size(); ++i ) {
      is_empty &= _q[i].empty();
    }
    return is_empty;
  }

  bool empty(const size_t i) {
    ASSERT(i < _q.size());
    return _q[i].empty();
  }

  void push(const T elem, const size_t i) {
    ASSERT(i < _q.size());
    return _q[i].push(elem);
  }

  T front(const size_t i) {
    ASSERT(i < _q.size());
    return _q[i].front();
  }

  void pop(const size_t i) {
    ASSERT(i < _q.size());
    return _q[i].pop();
  }

 private:
  vec<parallel::scalable_queue<T>> _q;
};
}

void ParallelConstruction::determineDistanceFromCut(const PartitionedHypergraph& phg,
                                                    const whfc::Node source,
                                                    const whfc::Node sink,
                                                    const PartitionID block_0,
                                                    const PartitionID block_1,
                                                    const vec<HypernodeID>& whfc_to_node) {
  _hfc.cs.borderNodes.distance.distance.assign(_flow_hg.numNodes(), whfc::HopDistance(0));
  _visited_hns.resize(_flow_hg.numNodes() + _flow_hg.numHyperedges());
  _visited_hns.reset();
  _visited_hns.set(source, true);
  _visited_hns.set(sink, true);

  // Initialize bfs queue with vertices contained in cut hyperedges
  size_t q_idx = 0;
  vec<BFSQueue<whfc::Node>> q(2, BFSQueue<whfc::Node>(_context.shared_memory.num_threads));
  tbb::parallel_for(0UL, _cut_hes.size(), [&](const size_t i) {
    const whfc::Hyperedge he = _flow_hg.originalHyperedgeID(_cut_hes[i].bucket, _cut_hes[i].e);
    for ( const whfc::FlowHypergraph::Pin& pin : _flow_hg.pinsOf(he) ) {
      if ( _visited_hns.compare_and_set_to_true(pin.pin) ) {
        q[q_idx].push(pin.pin, tbb::this_task_arena::current_thread_index());
      }
    }
    _visited_hns.set(_flow_hg.numNodes() + he, true);
  });

  // Perform BFS to determine distance of each vertex from cut
  std::vector<whfc::HopDistance>& distance = _hfc.cs.borderNodes.distance.distance;
  whfc::HopDistance dist(1);
  whfc::HopDistance max_dist_source(0);
  whfc::HopDistance max_dist_sink(0);
  while ( !q[q_idx].empty() ) {
    bool reached_source_side = false;
    bool reached_sink_side = false;
    tbb::parallel_for(0UL, _context.shared_memory.num_threads, [&](const size_t idx) {
      while ( !q[q_idx].empty(idx) ) {
        whfc::Node u = q[q_idx].front(idx);
        q[q_idx].pop(idx);
        const PartitionID block_of_u = phg.partID(whfc_to_node[u]);
        if ( block_of_u == block_0 ) {
          distance[u] = -dist;
          reached_source_side = true;
        } else if ( block_of_u == block_1 ) {
          distance[u] = dist;
          reached_sink_side = true;
        }

        for ( const whfc::FlowHypergraph::InHe& in_he : _flow_hg.hyperedgesOf(u) ) {
          const whfc::Hyperedge he = in_he.e;
          if ( _visited_hns.compare_and_set_to_true(_flow_hg.numNodes() + he) ) {
            for ( const whfc::FlowHypergraph::Pin& pin : _flow_hg.pinsOf(he) ) {
              if ( _visited_hns.compare_and_set_to_true(pin.pin) ) {
                q[1 - q_idx].push(pin.pin, idx);
              }
            }
          }
        }
      }
    });

    if ( reached_source_side ) max_dist_source = dist;
    if ( reached_sink_side ) max_dist_sink = dist;

    ASSERT(q[q_idx].empty());
    q_idx = 1 - q_idx;
    ++dist;
  }
  distance[source] = -(max_dist_source + 1);
  distance[sink] = max_dist_sink + 1;
}

} // namespace mt_kahypar
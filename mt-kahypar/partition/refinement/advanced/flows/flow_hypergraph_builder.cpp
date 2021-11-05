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

#include "mt-kahypar/partition/refinement/advanced/flows/flow_hypergraph_builder.h"

#include "tbb/blocked_range.h"
#include "tbb/parallel_invoke.h"
#include "tbb/parallel_reduce.h"
#include "tbb/parallel_for.h"

namespace mt_kahypar {

// ####################### Sequential Construction #######################

void FlowHypergraphBuilder::finalize() {
  if( !finishHyperedge() )	{ //finish last open hyperedge
    // maybe the last started hyperedge has zero/one pins and thus we still use the
    // previous sentinel. was never a bug, since that capacity is never read
    hyperedges.back().capacity = 0;
  }

  total_node_weight = whfc::NodeWeight(0);
  for (whfc::Node u : nodeIDs()) {
    nodes[u+1].first_out += nodes[u].first_out;
    total_node_weight += nodes[u].weight;
  }

  incident_hyperedges.resize(numPins());
  for (whfc::Hyperedge e : hyperedgeIDs()) {
    for (auto pin_it = beginIndexPins(e); pin_it != endIndexPins(e); pin_it++) {
      Pin& p = pins[pin_it];
      //destroy first_out temporarily and reset later
      whfc::InHeIndex ind_he = nodes[p.pin].first_out++;
      incident_hyperedges[ind_he] = { e, whfc::Flow(0), pin_it };
      //set iterator for incident hyperedge -> its position in incident_hyperedges of the node
      p.he_inc_iter = ind_he;
    }
  }

  for (whfc::Node u(numNodes()-1); u > 0; u--) {
    nodes[u].first_out = nodes[u-1].first_out;	//reset temporarily destroyed first_out
  }
  nodes[0].first_out = whfc::InHeIndex(0);

  _finalized = true;
}

bool FlowHypergraphBuilder::finishHyperedge() {
  if (currentHyperedgeSize() == 1) {
    removeLastPin();
  }

  if (currentHyperedgeSize() > 0) {
    pins_sending_flow.emplace_back(hyperedges.back().first_out, hyperedges.back().first_out);
    hyperedges.push_back({whfc::PinIndex::fromOtherValueType(numPins()), whfc::Flow(0), whfc::Flow(0)});//sentinel
    pins_receiving_flow.emplace_back(hyperedges.back().first_out, hyperedges.back().first_out);
    return true;
  }
  return false;
}

// ####################### Parallel Construction #######################

void FlowHypergraphBuilder::allocateHyperedgesAndPins(const size_t num_hyperedges,
                                                      const size_t num_pins) {
  tbb::parallel_invoke([&] {
    hyperedges.assign(num_hyperedges + 1, HyperedgeData {
      whfc::PinIndex::Invalid(), whfc::Flow(0), whfc::Flow(0) });
  }, [&] {
    pins.assign(num_pins, Pin { whfc::Node::Invalid(), whfc::InHeIndex::Invalid() });
  }, [&] {
    pins_sending_flow.assign(num_hyperedges, PinIndexRange {
      whfc::PinIndex::Invalid(), whfc::PinIndex::Invalid() });
  }, [&] {
    pins_receiving_flow.assign(num_hyperedges, PinIndexRange {
      whfc::PinIndex::Invalid(), whfc::PinIndex::Invalid() });
  });
}

void FlowHypergraphBuilder::finalizeHyperedges() {
  for ( size_t i = 1; i < _tmp_csr_buckets.size(); ++i ) {
    _tmp_csr_buckets[i]._global_start_he =
      _tmp_csr_buckets[i - 1]._global_start_he + _tmp_csr_buckets[i - 1]._num_hes;
    _tmp_csr_buckets[i]._global_start_pin_idx =
      _tmp_csr_buckets[i - 1]._global_start_pin_idx + _tmp_csr_buckets[i - 1]._num_pins;
  }

  tbb::parallel_for(0UL, _tmp_csr_buckets.size(), [&](const size_t idx) {
    _tmp_csr_buckets[idx].copyDataToFlowHypergraph(hyperedges, pins);
  });

  const size_t num_hyperedges =
    _tmp_csr_buckets.back()._global_start_he + _tmp_csr_buckets.back()._num_hes;
  const size_t num_pins =
    _tmp_csr_buckets.back()._global_start_pin_idx + _tmp_csr_buckets.back()._num_pins;
  resizeHyperedgesAndPins(num_hyperedges, num_pins);
  hyperedges.emplace_back( HyperedgeData { whfc::PinIndex(num_pins), whfc::Flow(0), whfc::Flow(0) } ); // sentinel
  tbb::parallel_for(0UL, num_hyperedges, [&](const size_t i) {
    pins_sending_flow[i].__begin = hyperedges[i].first_out;
    pins_sending_flow[i].__end = hyperedges[i].first_out;
    pins_receiving_flow[i].__begin = hyperedges[i + 1].first_out;
    pins_receiving_flow[i].__end = hyperedges[i + 1].first_out;
  });
}

void FlowHypergraphBuilder::finalizeParallel() {
  ASSERT(verifyParallelConstructedHypergraph(), "Parallel construction failed!");

  maxHyperedgeCapacity = tbb::parallel_reduce(
    tbb::blocked_range<size_t>(0UL, hyperedges.size()), whfc::Flow(0),
    [&](const tbb::blocked_range<size_t>& range, whfc::Flow init) {
      whfc::Flow max_capacity = init;
      for (size_t i = range.begin(); i < range.end(); ++i) {
        max_capacity = std::max(max_capacity, hyperedges[i].capacity);
      }
      return max_capacity;
    }, [](const whfc::Flow& lhs, const whfc::Flow& rhs) {
      return std::max(lhs, rhs);
    });

  total_node_weight = whfc::NodeWeight(0);
  for (whfc::Node u : nodeIDs()) {
    nodes[u+1].first_out += nodes[u].first_out;
    total_node_weight += nodes[u].weight;
  }

  incident_hyperedges.resize(numPins());
  for (whfc::Hyperedge e : hyperedgeIDs()) {
    for (auto pin_it = beginIndexPins(e); pin_it != endIndexPins(e); pin_it++) {
      Pin& p = pins[pin_it];
      //destroy first_out temporarily and reset later
      whfc::InHeIndex ind_he = nodes[p.pin].first_out++;
      incident_hyperedges[ind_he] = { e, whfc::Flow(0), pin_it };
      //set iterator for incident hyperedge -> its position in incident_hyperedges of the node
      p.he_inc_iter = ind_he;
    }
  }

  for (whfc::Node u(numNodes()-1); u > 0; u--) {
    nodes[u].first_out = nodes[u-1].first_out;	//reset temporarily destroyed first_out
  }
  nodes[0].first_out = whfc::InHeIndex(0);

  _finalized = true;
}


void FlowHypergraphBuilder::resizeHyperedgesAndPins(const size_t num_hyperedges,
                                                    const size_t num_pins) {
  ASSERT(num_hyperedges <= hyperedges.size());
  ASSERT(num_pins <= pins.size());
  hyperedges.resize(num_hyperedges);
  pins.resize(num_pins);
  pins_sending_flow.resize(num_hyperedges);
  pins_receiving_flow.resize(num_hyperedges);
}

bool FlowHypergraphBuilder::verifyParallelConstructedHypergraph() {
  size_t num_pins = 0;
  for ( size_t i = 0; i < numNodes(); ++i ) {
    if ( nodes[i].weight == 0 ) {
      LOG << "Node" << i << "has zero weight!";
      return false;
    }
    num_pins += nodes[i].first_out;
  }
  num_pins += nodes.back().first_out; // sentinel

  if ( num_pins != numPins() ) {
    LOG << "[Node Degrees] Expected number of pins =" << numPins() << ", Actual =" << num_pins;
    return false;
  }

  for ( size_t i = 0; i < pins.size(); ++i ) {
    if ( pins[i].pin == whfc::Node::Invalid() ) {
      LOG << "Pin at index" << i << "not assigned";
      return false;
    }
  }

  size_t previous_end = 0;
  num_pins = 0;
  for ( size_t i = 0; i < hyperedges.size() - 1; ++i ) {
    size_t current_start = hyperedges[i].first_out;
    size_t current_end = pins_receiving_flow[i].__end;
    if ( current_end - current_start <= 1 ) {
      LOG << "Hyperedge of size one contained";
      return false;
    }

    if ( current_start != previous_end ) {
      LOG << "Gap or intersection in hyperedge incidence array!";
      return false;
    }
    num_pins += ( current_end - current_start );
    previous_end = current_end;
  }

  if ( num_pins != numPins() ) {
    LOG << "[Edge Sizes] Expected number of pins =" << numPins() << ", Actual =" << num_pins;
    return false;
  }

  return true;
}

// ####################### Common Functions #######################

void FlowHypergraphBuilder::clear() {
  _finalized = false;
  _numPinsAtHyperedgeStart = 0;
  maxHyperedgeCapacity = 0;

  nodes.clear();
  hyperedges.clear();
  pins.clear();
  incident_hyperedges.clear();
  pins_sending_flow.clear();
  pins_receiving_flow.clear();
  total_node_weight = whfc::NodeWeight(0);
  sends_multiplier = 1;
  receives_multiplier = -1;

  //sentinels
  nodes.push_back({whfc::InHeIndex(0), whfc::NodeWeight(0)});
  hyperedges.push_back({whfc::PinIndex(0), whfc::Flow(0), whfc::Flow(0)});
}

}
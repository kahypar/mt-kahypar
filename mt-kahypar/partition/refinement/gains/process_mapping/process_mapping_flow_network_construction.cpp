/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/partition/refinement/gains/process_mapping/process_mapping_flow_network_construction.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/process_mapping/process_graph.h"

namespace mt_kahypar {

namespace {
HyperedgeWeight cap(const ProcessMappingCapacityAggregator aggregator,
                    const HyperedgeWeight edge_weight,
                    const HyperedgeWeight current_distance,
                    const HyperedgeWeight distance_without_block_0,
                    const HyperedgeWeight distance_without_block_1) {
  const HyperedgeWeight gain_0 = (current_distance - distance_without_block_0) * edge_weight;
  const HyperedgeWeight gain_1 = (current_distance - distance_without_block_1) * edge_weight;
  switch ( aggregator ) {
    case ProcessMappingCapacityAggregator::maximum: return std::max(gain_0, gain_1);
    case ProcessMappingCapacityAggregator::minimum: return std::min(gain_0, gain_1);
    case ProcessMappingCapacityAggregator::average: return (gain_0 + gain_1) / 2;
    case ProcessMappingCapacityAggregator::UNDEFINED:
      ERR("No valid capacity aggregator provided" << V(aggregator));
      return 0;
  }
  return 0;
}
} // namespace

template<typename PartitionedHypergraph>
HyperedgeWeight ProcessMappingFlowNetworkConstruction::capacity(const PartitionedHypergraph& phg,
                                                                const Context& context,
                                                                const HyperedgeID he,
                                                                const PartitionID block_0,
                                                                const PartitionID block_1)  {
  ASSERT(phg.hasProcessGraph());
  const ProcessGraph& process_graph = *phg.processGraph();
  const HyperedgeWeight edge_weight = phg.edgeWeight(he);
  const HypernodeID pin_count_block_0 = phg.pinCountInPart(he, block_0);
  const HypernodeID pin_count_block_1 = phg.pinCountInPart(he, block_1);
  ds::Bitset& connectivity_set = phg.deepCopyOfConnectivitySet(he);
  if ( pin_count_block_0 > 0 && pin_count_block_1 == 0 ) {
    // Hyperedge is non-cut in the flow network.
    // Thus, we set its capacity to gain if we would make it a cut edge (adding block 1)
    connectivity_set.set(block_0);
    return (process_graph.distanceWithBlock(connectivity_set, block_1) -
      process_graph.distance(connectivity_set)) * edge_weight;
  } else if ( pin_count_block_0 == 0 && pin_count_block_1 > 0 ) {
    // Hyperedge is non-cut in the flow network.
    // Thus, we set its capacity to gain if we would make it a cut edge (adding block 0)
    connectivity_set.set(block_1);
    return (process_graph.distanceWithBlock(connectivity_set, block_0) -
      process_graph.distance(connectivity_set)) * edge_weight;
  } else {
    // Hyperedge is cut in the flow network.
    // The gain for making the hyperedge non-cut depends on whether we
    // remove block 0 or 1 from the connectivity set. Thus, we compute both gains
    // and set the capacity of the edge to the result of an aggregation operation
    // over both gains (minimum, maximum or average)
    connectivity_set.set(block_0);
    connectivity_set.set(block_1);
    const HyperedgeWeight current_distance = process_graph.distance(connectivity_set);
    const HyperedgeWeight distance_without_block_0 = process_graph.distanceWithoutBlock(connectivity_set, block_0);
    const HyperedgeWeight distance_without_block_1 = process_graph.distanceWithoutBlock(connectivity_set, block_1);
    return cap(context.refinement.flows.capacity_aggregator, edge_weight,
      current_distance, distance_without_block_0, distance_without_block_1);
  }
}

namespace {
#define PROCESS_MAPPING_CAPACITY(X) HyperedgeWeight ProcessMappingFlowNetworkConstruction::capacity(  \
  const X&, const Context&, const HyperedgeID, const PartitionID, const PartitionID)
}

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(PROCESS_MAPPING_CAPACITY)

}  // namespace mt_kahypar

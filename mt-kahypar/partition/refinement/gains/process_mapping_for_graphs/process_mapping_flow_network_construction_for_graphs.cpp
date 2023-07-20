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

#include "mt-kahypar/partition/refinement/gains/process_mapping_for_graphs/process_mapping_flow_network_construction_for_graphs.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/process_mapping/process_graph.h"

namespace mt_kahypar {

template<typename PartitionedHypergraph>
HyperedgeWeight GraphProcessMappingFlowNetworkConstruction::capacity(const PartitionedHypergraph& phg,
                                                                     const Context&,
                                                                     const HyperedgeID he,
                                                                     const PartitionID block_0,
                                                                     const PartitionID block_1)  {
  ASSERT(phg.hasProcessGraph());
  const ProcessGraph& process_graph = *phg.processGraph();
  const HyperedgeWeight edge_weight = phg.edgeWeight(he);
  const HypernodeID u = phg.edgeSource(he);
  const HypernodeID v = phg.edgeTarget(he);
  const PartitionID block_of_u = phg.partID(u);
  const PartitionID block_of_v = phg.partID(v);

  if ( ( block_of_u == block_0 || block_of_u == block_1 ) &&
       ( block_of_v == block_0 || block_of_v == block_1 )  ) {
    // Both endpoints of the edge are either contained in block 0 or 1.
    // Removing the edge from the cut or making it a cut edge has the
    // following gain:
    return process_graph.distance(block_0, block_1) * edge_weight;
  } else {
    // In this case, only one node is contained in the flow problem and the other
    // node is part of another block different from block_0 and block_1.
    // Here, we set the capacity to difference in the steiner tree metric
    // if we would replace block_0 with block_1.
    PartitionID other_block = kInvalidPartition;
    if ( block_of_u == block_0 || block_of_v == block_0 ) {
      other_block = block_of_u == block_0 ? block_of_v : block_of_u;
    } else if ( block_of_u == block_1 || block_of_v == block_1 ) {
      other_block = block_of_u == block_1 ? block_of_v : block_of_u;
    } else {
      // Can happen due to concurrent node moves applied by other flow problems
      return 0;
    }
    ASSERT(other_block != kInvalidPartition);
    const HyperedgeWeight current_distance = process_graph.distance(block_0, other_block);
    const HyperedgeWeight distance_with_block_1 = process_graph.distance(block_1, other_block);
    return std::abs(current_distance - distance_with_block_1) * edge_weight;
  }
  return 0;
}

template<typename PartitionedHypergraph>
bool GraphProcessMappingFlowNetworkConstruction::connectToSource(const PartitionedHypergraph& partitioned_hg,
                                                                 const HyperedgeID he,
                                                                 const PartitionID block_0,
                                                                 const PartitionID block_1) {
  ASSERT(partitioned_hg.hasProcessGraph());
  const ProcessGraph& process_graph = *partitioned_hg.processGraph();
  const HypernodeID u = partitioned_hg.edgeSource(he);
  const HypernodeID v = partitioned_hg.edgeTarget(he);
  const PartitionID block_of_u = partitioned_hg.partID(u);
  const PartitionID block_of_v = partitioned_hg.partID(v);
  if ( block_of_u == block_0 || block_of_v == block_0 ) {
    PartitionID other_block = block_of_u == block_0 ? block_of_v : block_of_u;
    const HyperedgeWeight current_distance = process_graph.distance(block_0, other_block);
    const HyperedgeWeight distance_block_1 = process_graph.distance(block_1, other_block);
    if ( other_block != block_0 && other_block != block_1 && current_distance < distance_block_1 ) {
      // Moving the node from block_0 to block_1 would worsen the steiner tree metric,
      // even though the edge is still cut afterwards. To model this percurlarity in the flow network,
      // we add the corresponding edge to the source.
      return true;
    }
  }
  if ( block_of_u == block_1 || block_of_v == block_1 ) {
    PartitionID other_block = block_of_u == block_1 ? block_of_v : block_of_u;
    const HyperedgeWeight current_distance = process_graph.distance(block_1, other_block);
    const HyperedgeWeight distance_block_0 = process_graph.distance(block_0, other_block);
    if ( other_block != block_0 && other_block != block_1 && current_distance > distance_block_0 ) {
      // Moving the node from block_1 to block_0 would improve the steiner tree metric,
      // even though the edge is still cut afterwards. To model this percurlarity in the flow network,
      // we add the corresponding edge to the source.
      return true;
    }
  }
  return false;
}


template<typename PartitionedHypergraph>
bool GraphProcessMappingFlowNetworkConstruction::connectToSink(const PartitionedHypergraph& partitioned_hg,
                                                               const HyperedgeID he,
                                                               const PartitionID block_0,
                                                               const PartitionID block_1) {
  ASSERT(partitioned_hg.hasProcessGraph());
  const ProcessGraph& process_graph = *partitioned_hg.processGraph();
  const HypernodeID u = partitioned_hg.edgeSource(he);
  const HypernodeID v = partitioned_hg.edgeTarget(he);
  const PartitionID block_of_u = partitioned_hg.partID(u);
  const PartitionID block_of_v = partitioned_hg.partID(v);
  if ( block_of_u == block_1 || block_of_v == block_1 ) {
    PartitionID other_block = block_of_u == block_1 ? block_of_v : block_of_u;
    const HyperedgeWeight current_distance = process_graph.distance(block_1, other_block);
    const HyperedgeWeight distance_block_1 = process_graph.distance(block_0, other_block);
    if ( other_block != block_0 && other_block != block_1 && current_distance < distance_block_1 ) {
      // Moving the node from block_1 to block_0 would worsen the steiner tree metric,
      // even though the edge is still cut afterwards. To model this percurlarity in the flow network,
      // we add the corresponding edge to the sink.
      return true;
    }
  }
  if ( block_of_u == block_0 || block_of_v == block_0 ) {
    PartitionID other_block = block_of_u == block_0 ? block_of_v : block_of_u;
    const HyperedgeWeight current_distance = process_graph.distance(block_0, other_block);
    const HyperedgeWeight distance_block_1 = process_graph.distance(block_1, other_block);
    if ( other_block != block_0 && other_block != block_1 && current_distance > distance_block_1 ) {
      // Moving the node from block_0 to block_1 would improve the steiner tree metric,
      // even though the edge is still cut afterwards. To model this percurlarity in the flow network,
      // we add the corresponding edge to the sink.
      return true;
    }
  }
  return false;
}

template<typename PartitionedHypergraph>
bool GraphProcessMappingFlowNetworkConstruction::isCut(const PartitionedHypergraph& partitioned_hg,
                                                       const HyperedgeID he,
                                                       const PartitionID block_0,
                                                       const PartitionID block_1) {
  ASSERT(partitioned_hg.hasProcessGraph());
  const ProcessGraph& process_graph = *partitioned_hg.processGraph();
  const HypernodeID u = partitioned_hg.edgeSource(he);
  const HypernodeID v = partitioned_hg.edgeTarget(he);
  const PartitionID block_of_u = partitioned_hg.partID(u);
  const PartitionID block_of_v = partitioned_hg.partID(v);
  if ( block_of_u == block_1 || block_of_v == block_1 ) {
    PartitionID other_block = block_of_u == block_1 ? block_of_v : block_of_u;
    const HyperedgeWeight current_distance = process_graph.distance(block_1, other_block);
    const HyperedgeWeight distance_block_0 = process_graph.distance(block_0, other_block);
    if ( other_block != block_0 && other_block != block_1 && current_distance > distance_block_0 ) {
      // Moving the node contained in the flow problem to the other block would improve the
      // steiner tree metric, even though the edge would be still cut.
      // Thus, we consider it as a cut edge.
      return true;
    }
  }
  if ( block_of_u == block_0 || block_of_v == block_0 ) {
    PartitionID other_block = block_of_u == block_0 ? block_of_v : block_of_u;
    const HyperedgeWeight current_distance = process_graph.distance(block_0, other_block);
    const HyperedgeWeight distance_block_1 = process_graph.distance(block_1, other_block);
    if ( other_block != block_0 && other_block != block_1 && current_distance > distance_block_1 ) {
      // Same as the previous case
      return true;
    }
  }
  return false;
}

namespace {
#define PROCESS_MAPPING_CAPACITY(X) HyperedgeWeight GraphProcessMappingFlowNetworkConstruction::capacity(  \
  const X&, const Context&, const HyperedgeID, const PartitionID, const PartitionID)
#define PROCESS_MAPPING_CONNECT_TO_SOURCE(X) bool GraphProcessMappingFlowNetworkConstruction::connectToSource(  \
  const X&, const HyperedgeID, const PartitionID, const PartitionID)
#define PROCESS_MAPPING_CONNECT_TO_SINK(X) bool GraphProcessMappingFlowNetworkConstruction::connectToSink(  \
  const X&, const HyperedgeID, const PartitionID, const PartitionID)
#define PROCESS_MAPPING_IS_CUT(X) bool GraphProcessMappingFlowNetworkConstruction::isCut(  \
  const X&, const HyperedgeID, const PartitionID, const PartitionID)
}

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(PROCESS_MAPPING_CAPACITY)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(PROCESS_MAPPING_CONNECT_TO_SOURCE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(PROCESS_MAPPING_CONNECT_TO_SINK)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(PROCESS_MAPPING_IS_CUT)

}  // namespace mt_kahypar

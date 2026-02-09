#pragma once

#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/partition/refinement/i_rebalancer.h"
#include "mt-kahypar/datastructures/priority_queue.h"
#include "mt-kahypar/utils/cast.h"

namespace mt_kahypar::constraints {

using Key = std::tuple<HypernodeID, HypernodeID, HypernodeID>;

template<typename Comparator = std::less<Key>, uint32_t arity = 4>
using PQ = ds::Heap<Key, HypernodeID, Comparator, arity>;

template<typename PartitionedHypergraph>
bool verifyConstraints(const PartitionedHypergraph& partitioned_hg, const Context& context) {
  vec<std::pair<HypernodeID, HypernodeID>> constraints;
  io::readNegativeConstraintsFile(context.partition.negative_constraints_filename, constraints);
  for (std::pair<HypernodeID, HypernodeID> constraint : constraints) {
    if (partitioned_hg.partID(constraint.first) == partitioned_hg.partID(constraint.second)) {
      return false;
    }
  }
  return true;
}

template<typename PartitionedHypergraph>
HypernodeID incientNodesInSamePart(const PartitionedHypergraph& partitioned_hg,
                                        const HypernodeID& node_id) {
  HypernodeID node;
  HypernodeID num_nodes = 0;
  if(partitioned_hg.fixedVertexSupport().getConstraintIdFromHypergraphId(node_id, node) &&
    partitioned_hg.partID(node_id) != kInvalidPartition ) {

    PartitionID part = partitioned_hg.partID(node_id);
    const ds::DynamicGraph& constraint_graph = partitioned_hg.fixedVertexSupport().getConstraintGraph();
    for (HypernodeID incident_node : constraint_graph.incidentNodes(node)) {
      PartitionID pid = partitioned_hg.partID(constraint_graph.nodeWeight(incident_node));
      if (pid == part) num_nodes++;
    }
  }
  return num_nodes;
}

template<typename PartitionedHypergraph>
PartitionID getBestCutPartition(const bool must_cut_be_positive,
                                      const HypernodeID& node_id,
                                      const PartitionID& current_partition,
                                      const vec<bool>& is_partition_invalid,
                                      gain_cache_t& gain_cache) {
  PartitionID best_partition = current_partition;
  const PartitionID num_partitons = is_partition_invalid.size();
  HyperedgeWeight max_gain = std::numeric_limits<HyperedgeWeight>::min();

  for (PartitionID partition = 0; partition < num_partitons; partition++) {
    if (is_partition_invalid[partition]) continue;
    HyperedgeWeight gain = GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>(
      [&](const auto& cache){
        return cache.gain(node_id, current_partition, partition);
      },
      gain_cache
    );
    if (must_cut_be_positive && gain < 0) continue;
    if (gain > max_gain) {
      max_gain = gain;
      best_partition = partition;
    }
  }
  return best_partition;
}

template<typename PartitionedHypergraph>
void descendingConstraintDegree(PartitionedHypergraph& partitioned_hg,
                                const Context& context,
                                gain_cache_t& gain_cache) {
  unused(context);
  GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>(
    [&](auto& cache){
        cache.initializeGainCache(partitioned_hg);
    },
    gain_cache
  );
  auto delta_func =
  [&partitioned_hg, &gain_cache](const SynchronizedEdgeUpdate& sync_update) {
    GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>(
      [&](auto& cache) {
        cache.deltaGainUpdate(partitioned_hg, sync_update);
      },
      gain_cache
    );
  };
  // initialize PQ
  const ds::DynamicGraph& constraint_graph = partitioned_hg.fixedVertexSupport().getConstraintGraph();
  std::vector<PosT> positions(constraint_graph.numNodes(), invalid_position);
  PQ<std::less<Key>, 4> heap(positions.data(), positions.size());
  for (auto node : constraint_graph.nodes()) {
    heap.insert(node, {incientNodesInSamePart(partitioned_hg, node), constraint_graph.nodeDegree(node), node});
  }
  
  while(!heap.empty()){
    HypernodeID node = heap.top();
    heap.deleteTop();
    HypernodeID node_id = HypernodeID(constraint_graph.nodeWeight(node));
    PartitionID partition_id = partitioned_hg.partID(node_id);
    vec<bool> invalid_partitions(partitioned_hg.k(), false);

    for (HypernodeID incident_node : constraint_graph.incidentNodes(node)) {
        PartitionID incident_partition_id = partitioned_hg.partID(constraint_graph.nodeWeight(incident_node));
        invalid_partitions[incident_partition_id] = true;
    }
    PartitionID new_partition_id = getBestCutPartition<PartitionedHypergraph>(false, node_id, partition_id, invalid_partitions, gain_cache);
    if (new_partition_id != partition_id) {
      partitioned_hg.changeNodePart(node_id, partition_id, new_partition_id, delta_func);
    }
  }
}

template<typename PartitionedHypergraph>
void postprocessNegativeConstraints(PartitionedHypergraph& partitioned_hg,
                                    const Context& context) {
  gain_cache_t gain_cache = GainCachePtr::constructGainCache(context);
  std::unique_ptr<IRebalancer> rebalancer = RebalancerFactory::getInstance().createObject(
      context.refinement.rebalancing.algorithm, partitioned_hg.initialNumNodes(), context, gain_cache);
  descendingConstraintDegree(partitioned_hg, context, gain_cache);

  Metrics metrics { metrics::quality(partitioned_hg, context), metrics::imbalance(partitioned_hg, context) };
  mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(partitioned_hg);
  rebalancer->initialize(phg);
  rebalancer->refine(phg, {}, metrics, 0.0);
  LOG << "-------------- stats after postprocessing --------------";
  LOG << "km1       ="<< metrics::quality(partitioned_hg, context);
  LOG << "Imbalance ="<<metrics::imbalance(partitioned_hg, context);
  LOG << (verifyConstraints(partitioned_hg, context)? "Constrains were respected from balancer" : "!!! Balancer destroyed constrains !!!");
  LOG << "";
  GainCachePtr::deleteGainCache(gain_cache);
}

} // namespace mt_kahypar::constraints
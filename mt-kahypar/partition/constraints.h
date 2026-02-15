#pragma once

#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/partition/refinement/i_rebalancer.h"
#include "mt-kahypar/utils/cast.h"

namespace mt_kahypar::constraints {

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
PartitionID getBestCutPartition(const HypernodeID& node_id,
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
    if (gain > max_gain) {
      max_gain = gain;
      best_partition = partition;
    }
  }
  return best_partition;
}

template<typename PartitionedHypergraph>
void frontToBackConstraints(PartitionedHypergraph& partitioned_hg,
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

  /**
   * For each node in constraint graph 
   * -> get neighbors
   * -> get node_ids in the PartitionedHypergraph from the node_weights in the constraint graph
   * -> get all partitionIDs 
   * -> move node_id in different partition if nessesary
   */
  const ds::DynamicGraph& constraint_graph = partitioned_hg.fixedVertexSupport().getConstraintGraph();

  for ( const auto& node : constraint_graph.nodes()) {
      const HypernodeID node_id = HypernodeID(constraint_graph.nodeWeight(node));
      const PartitionID partition_id = partitioned_hg.partID(node_id);
      vec<bool> invalid_partitions(partitioned_hg.k(), false);

      for (HypernodeID incident_node : constraint_graph.incidentNodes(node)) {
          PartitionID incident_partition_id = partitioned_hg.partID(constraint_graph.nodeWeight(incident_node));
          invalid_partitions[incident_partition_id] = true;
      }
      PartitionID new_partition_id = getBestCutPartition<PartitionedHypergraph>(node_id, partition_id, invalid_partitions, gain_cache);
      if ( invalid_partitions[partition_id] && new_partition_id != partition_id ) {
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

  frontToBackConstraints(partitioned_hg, context, gain_cache);

  LOG << "";
  LOG << "Verify if constraints are respected:";
  LOG << "";
  LOG << (verifyConstraints(partitioned_hg, context)? "Constrains were respected from partitioner" : "!!! Partitioner destroyed constrains !!!");

  Metrics metrics { metrics::quality(partitioned_hg, context), metrics::imbalance(partitioned_hg, context) };
  mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(partitioned_hg);
  rebalancer->initialize(phg);
  rebalancer->refine(phg, {}, metrics, 0.0);
  GainCachePtr::deleteGainCache(gain_cache);
  LOG << (verifyConstraints(partitioned_hg, context)? "Constrains were respected from balancer" : "!!! Balancer destroyed constrains !!!");
  LOG << "";
}

} // namespace mt_kahypar::constraints
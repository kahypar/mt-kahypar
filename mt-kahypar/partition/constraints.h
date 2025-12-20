#pragma once

#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/fixed_vertex_support.h"
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

template<typename Hypergraph>
bool verifyConstraints(const Hypergraph& hg) {
  const ds::FixedVertexSupport<Hypergraph>& fixed_vertex_support = hg.fixedVertexSupport();
  for (const auto& constraint : fixed_vertex_support.getConstraints()) {
    HypernodeID u;
    HypernodeID v;
    if (fixed_vertex_support.getConstraintIdFromHypergraphId(constraint.first, u) && fixed_vertex_support.getConstraintIdFromHypergraphId(constraint.second, v)) {
      // constraint between u and v exists, but they are contracted together
      if (u == v) return false;
    } else {
      // ohne of the nodes was not found in
      throw std::logic_error("Node is in constraints but not in _hypergraph_id_to_graph_id mapping");
      return false;}
  }
  return true;
}

template<typename PartitionedHypergraph>
  PartitionID getLowestWeightPartition(const PartitionedHypergraph& partitioned_hg,
                                  const Context& context,
                                  const HypernodeID& node_id,
                                  const vec<bool>& is_partition_invalid,
                                  const Km1GainCache& concrete_gain_cache) {
    PartitionID best_partition = partitioned_hg.partID(999);
    assert(is_partition_invalid[best_partition]);
    const PartitionID num_partitons = is_partition_invalid.size();
    const HypernodeWeight node_weight = partitioned_hg.nodeWeight(node_id);
    HyperedgeWeight max_gain = std::numeric_limits<HyperedgeWeight>::min();
    HypernodeWeight min_weight = std::numeric_limits<HypernodeWeight>::max();

    PartitionID fallback_partition = best_partition;
    HypernodeWeight fallback_weight = std::numeric_limits<HypernodeWeight>::max();

    for (PartitionID partition = 0; partition < num_partitons; partition++) {
      if (is_partition_invalid[partition]) continue;

      HypernodeWeight max_weight = context.partition.max_part_weights[partition];
      HypernodeWeight weight = partitioned_hg.partWeight(partition);
      HyperedgeWeight gain = concrete_gain_cache.benefitTerm(node_id, partition);
      if (weight < fallback_weight) {
        fallback_weight = weight;
        fallback_partition = partition;
      }

      if ((gain > max_gain || (gain == max_gain && weight < min_weight)) &&
          weight + node_weight <= max_weight) {
          max_gain = gain;
          min_weight = weight;
          best_partition = partition;
      }
    }
    if (max_gain == std::numeric_limits<HyperedgeWeight>::min()) {
      return fallback_partition;
    }
    return best_partition;
  }

template<typename PartitionedHypergraph>
void postprocessNegativeConstraints(PartitionedHypergraph& partitioned_hg,
                                    const Context& context) {
    /**
     * For each node in constraint graph 
     * -> get neighbors
     * -> get node_ids in the PartitionedHypergraph from the node_weights in the constraint graph
     * -> get all partitionIDs 
     * -> move node_id in different partition if nessesary
     */
    gain_cache_t gain_cache = GainCachePtr::constructGainCache(context);
    Km1GainCache& concrete_gain_cache = GainCachePtr::cast<Km1GainCache>(gain_cache);
    concrete_gain_cache.initializeGainCache(partitioned_hg);
    std::unique_ptr<IRebalancer> rebalancer = RebalancerFactory::getInstance().createObject(
        context.refinement.rebalancing.algorithm, partitioned_hg.initialNumNodes(), context, gain_cache);

    const ds::DynamicGraph& constraint_graph = partitioned_hg.fixedVertexSupport().getConstraintGraph();

    for ( const auto& node : constraint_graph.nodes()) {
        const HypernodeID node_id = HypernodeID(constraint_graph.nodeWeight(node));
        const PartitionID partition_id = partitioned_hg.partID(node_id);
        vec<bool> invalid_partitions(partitioned_hg.k(), false);

        for (HypernodeID incident_node : constraint_graph.incidentNodes(node)) {
            PartitionID incident_partition_id = partitioned_hg.partID(constraint_graph.nodeWeight(incident_node));
            invalid_partitions[incident_partition_id] = true;
        }
        if (invalid_partitions[partition_id]) {
            PartitionID new_partition_id = getLowestWeightPartition(partitioned_hg, context, node_id, invalid_partitions, concrete_gain_cache);
            partitioned_hg.changeNodePart(node_id,
                                        partition_id,
                                        new_partition_id);
        }
        // LOG << "Node nr: " << node_id;
        // LOG << "Partition id: " << partition_id;
        // LOG << (invalid_partitions[partition_id]? ("Moved to Partition:") : ("Stayed in:")) << partitioned_hg.partID(node_id);
        // LOG << "";
    }

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
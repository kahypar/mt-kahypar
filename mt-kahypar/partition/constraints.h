#pragma once

#include <unordered_set>
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/fixed_vertex_support.h"
#include "mt-kahypar/datastructures/priority_queue.h"
#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/partition/refinement/i_rebalancer.h"
#include "mt-kahypar/utils/cast.h"

namespace mt_kahypar::constraints {

// using Key = std::tuple<PartitionID, PartitionID, HypernodeID>;

// template<typename Comparator = std::less<Key>, uint32_t arity = 4>
// using PQ = ds::Heap<Key, HypernodeID, Comparator, arity>;

template<typename PartitionedHypergraph>
PartitionID constraintDegree(const PartitionedHypergraph& partitioned_hg,
                                        const HypernodeID& node_id) {
  HypernodeID node;
  PartitionID num_constraints = 0;
  if(partitioned_hg.fixedVertexSupport().getConstraintIdFromHypergraphId(node_id, node)) {
    vec<bool> is_partition_allowed(partitioned_hg.k(), true);
    const ds::DynamicGraph& constraint_graph = partitioned_hg.fixedVertexSupport().getConstraintGraph();
    for (HypernodeID incident_node : constraint_graph.incidentNodes(node)) {
      PartitionID pid = partitioned_hg.partID(constraint_graph.nodeWeight(incident_node));
      if (pid < 0 || pid >= partitioned_hg.k()) continue; // skip invalid part ids, not all nodes are already in part
      is_partition_allowed[pid] = false;
    }
    for (bool allowed : is_partition_allowed) {
      if (!allowed) num_constraints++;
    }
  }
  return num_constraints;
}

template<typename PartitionedHypergraph>
PartitionID numberOfDistinctNeighbors(const ds::DynamicGraph& constraint_graph, HypernodeID u) {
  PartitionID count = 0;
  std::unordered_set<HypernodeID> set;
  for (const auto& node : constraint_graph.incidentNodes(u)) {
    if (set.find(node) != set.end()) continue;
    else {
      set.insert(node);
      count++;
    }
  }
  return count;
}

template<typename PartitionedHypergraph>
bool constraintsMet(const PartitionedHypergraph& partitioned_hg, const Context& context) {
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
HypernodeID numBrokenConstraints(const PartitionedHypergraph& hg) {
  HypernodeID count = 0;
  for (const auto& constraint : hg.fixedVertexSupport().getConstraints()) {
    if (hg.partID(constraint.first) == hg.partID(constraint.second) && hg.partID(constraint.first) != kInvalidPartition) {
      count++;
    }
  }
  return count;
}

template<typename PartitionedHypergraph>
bool constraintsMet(const PartitionedHypergraph& hg) {
  return numBrokenConstraints(hg) == 0;
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
bool isNodeAllowedInPartition(const PartitionedHypergraph& partitioned_hg,
                            const HypernodeID& node_id,
                            const PartitionID part_id) {
  HypernodeID node;
  if(partitioned_hg.fixedVertexSupport().getConstraintIdFromHypergraphId(node_id, node)) {
    const ds::DynamicGraph& constraint_graph = partitioned_hg.fixedVertexSupport().getConstraintGraph();
    vec<HypernodeID> count_constraints_in_partition(partitioned_hg.k(), 0);
    for (HypernodeID incident_node : constraint_graph.incidentNodes(node)) {
      PartitionID pid = partitioned_hg.partID(constraint_graph.nodeWeight(incident_node));
      if (pid >= 0 && pid < partitioned_hg.k()){
        count_constraints_in_partition[pid]++;
      }
    }
    HypernodeID min_value = *std::min_element(count_constraints_in_partition.begin(), count_constraints_in_partition.end());
    if (count_constraints_in_partition[part_id] == min_value) return true;
    else return false;
  }
  return true;
}

template<typename PartitionedHypergraph>
bool isNodeAllowedInPartition(const PartitionedHypergraph& partitioned_hg,
                            const Context& context,
                            const HypernodeID& node_id,
                            const PartitionID part_id) {
  HypernodeID node;
  if(partitioned_hg.fixedVertexSupport().getConstraintIdFromHypergraphId(node_id, node)) {
    const ds::DynamicGraph& constraint_graph = partitioned_hg.fixedVertexSupport().getConstraintGraph();
    PartitionID upper_level_allowed_constraints = context.coarsening.allowed_constraint_degree;

    // only allow that amount of constraints in block k that can be later be put in different blocks
    HypernodeID allowed_constraints = upper_level_allowed_constraints / 2 - 1;
    if(part_id == 0 && upper_level_allowed_constraints % 2 == 1) {
      allowed_constraints++; // if there is an odd number of blocks left, block 0 gets one time more devided than block 1
    }

    HypernodeID constraint_count = 0;
    for (HypernodeID incident_node : constraint_graph.incidentNodes(node)) {
      const HypernodeID incident_node_id = constraint_graph.nodeWeight(incident_node);
      const PartitionID pid = partitioned_hg.partID(incident_node_id);
      if (pid == part_id) {
        constraint_count++;
        if (incientNodesInSamePart(partitioned_hg, incident_node_id) >= allowed_constraints) {
          //if current part is not allowed, choose part with least constraints
          return false;//isNodeAllowedInPartition(partitioned_hg, node_id, part_id);
        }
      }
    }
    return constraint_count <= allowed_constraints;
  }
  return true;
}

template<typename PartitionedHypergraph>
PartitionID allNodesAllowedNumberOfNeighbors(const PartitionedHypergraph& partitioned_hg) {
  // contraint nodes that get contracted must have less than k neighbors! otherwise they break the invariant and initial partitioning is impossible
  
  bool pass = true;
  PartitionID k = partitioned_hg.fixedVertexSupport().numBlocks();
  const ds::DynamicGraph& constraint_graph = partitioned_hg.fixedVertexSupport().getConstraintGraph();
  for(const auto& u : constraint_graph.nodes()){
    PartitionID count = numberOfDistinctNeighbors<PartitionedHypergraph>(constraint_graph, u);
    if(count >= k) {
      //LOG << "node" << constraint_graph.nodeWeight(u) << "has too much neighbors!" << count;
      pass = false;
      }
  }
  return pass;
}

template<typename PartitionedHypergraph>
PartitionID allNodesAllowedNumberOfNeighbors(const PartitionedHypergraph& partitioned_hg, PartitionID k) {
  // contraint nodes that get contracted must have less than k neighbors! otherwise they break the invariant and initial partitioning is impossible
  
  bool pass = true;
  const ds::DynamicGraph& constraint_graph = partitioned_hg.fixedVertexSupport().getConstraintGraph();
  for(const auto& u : constraint_graph.nodes()){
    PartitionID count = numberOfDistinctNeighbors<PartitionedHypergraph>(constraint_graph, u);
    if(count >= k) {
      LOG << "node" << constraint_graph.nodeWeight(u) << "has too much neighbors!" << count;
      pass = false;
      }
  }
  return pass;
}

template<typename PartitionedHypergraph>
PartitionID numberOfDistinctNeighbors(const ds::DynamicGraph& constraint_graph,
                            const HypernodeID& u,
                            const HypernodeID& v) {
  // contraint nodes that get contracted must have less than k neighbors! otherwise they break the invariant and initial partitioning is impossible
  PartitionID count = 0;
  std::unordered_set<HypernodeID> set;
  for (const auto& node : constraint_graph.incidentNodes(u)) {
    if (set.find(node) != set.end()) continue; // TODO: there are many double neighbors, how to remove them?
    else {
      set.insert(node);
      count++;
    }
  }
  for (const auto& node : constraint_graph.incidentNodes(v)) {
    if (set.find(node) != set.end()) continue;
    else {
      set.insert(node);
      count++;
    }
  }
  return count;
}



template<typename PartitionedHypergraph>
PartitionID getLowestWeightPartition(const PartitionedHypergraph& partitioned_hg,
                                      const Context& context,
                                      const HypernodeID& node_id,
                                      const PartitionID& current_partition,
                                      const vec<bool>& is_partition_invalid,
                                      const Km1GainCache& concrete_gain_cache) {
  PartitionID best_partition = kInvalidPartition;
  assert(is_partition_invalid[partitioned_hg.partID(node_id)]);
  const PartitionID num_partitons = is_partition_invalid.size();
  const HypernodeWeight node_weight = partitioned_hg.nodeWeight(node_id);
  HyperedgeWeight max_gain = std::numeric_limits<HyperedgeWeight>::min();
  HypernodeWeight min_weight = std::numeric_limits<HypernodeWeight>::max();

  PartitionID fallback_partition = kInvalidPartition;
  HypernodeWeight fallback_weight = std::numeric_limits<HypernodeWeight>::max();

  for (PartitionID partition = 0; partition < num_partitons; partition++) {
    if (is_partition_invalid[partition]) continue;

    HypernodeWeight max_weight = context.partition.max_part_weights[partition];
    HypernodeWeight weight = partitioned_hg.partWeight(partition);
    HyperedgeWeight gain = concrete_gain_cache.gain(node_id, current_partition, partition);
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

  for (int i = 0; i < 5; i++) {
    for ( const auto& node : constraint_graph.nodes()) {
        const HypernodeID node_id = HypernodeID(constraint_graph.nodeWeight(node));
        const PartitionID partition_id = partitioned_hg.partID(node_id);
        vec<bool> invalid_partitions(partitioned_hg.k(), false);

        for (HypernodeID incident_node : constraint_graph.incidentNodes(node)) {
            PartitionID incident_partition_id = partitioned_hg.partID(constraint_graph.nodeWeight(incident_node));
            invalid_partitions[incident_partition_id] = true;
        }
        // in first round just move with gain >= 0
        PartitionID new_partition_id = getBestCutPartition<PartitionedHypergraph>(false, node_id, partition_id, invalid_partitions, gain_cache);
        if ( new_partition_id != partition_id) {
          partitioned_hg.changeNodePart(node_id, partition_id, new_partition_id, delta_func);
        }
    }
  }
}

// template<typename PartitionedHypergraph>
// void descendingConstraintDegree(PartitionedHypergraph& partitioned_hg,
//                                 const Context& context,
//                                 gain_cache_t& gain_cache) {
//   unused(context);
//   GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>(
//     [&](auto& cache){
//         cache.initializeGainCache(partitioned_hg);
//     },
//     gain_cache
//   );
//   // initialize PQ
//   const ds::DynamicGraph& constraint_graph = partitioned_hg.fixedVertexSupport().getConstraintGraph();
//   std::vector<PosT> positions(constraint_graph.numNodes(), invalid_position);
//   PQ heap(positions.data(), positions.size());
//   for (auto node : constraint_graph.nodes()) {
//     HypernodeID node_id = constraint_graph.nodeWeight(node);
//     heap.insert(node, {constraintDegree(partitioned_hg, node), partitioned_hg.nodeDegree(node_id), node});
//   }
  
//   while(!heap.empty()){
//     HypernodeID node = heap.top();
//     heap.deleteTop();
//     HypernodeID node_id = HypernodeID(constraint_graph.nodeWeight(node));
//     PartitionID partition_id = partitioned_hg.partID(node_id);
//     vec<bool> invalid_partitions(partitioned_hg.k(), false);

//     for (HypernodeID incident_node : constraint_graph.incidentNodes(node)) {
//         PartitionID incident_partition_id = partitioned_hg.partID(constraint_graph.nodeWeight(incident_node));
//         invalid_partitions[incident_partition_id] = true;
//     }
//     PartitionID new_partition_id = getBestCutPartition<PartitionedHypergraph>(false, node_id, partition_id, invalid_partitions, gain_cache);
//     if (new_partition_id != partition_id) {
//       partitioned_hg.changeNodePart(node_id, partition_id, new_partition_id);
//     }
//   }
// }

template<typename PartitionedHypergraph>
void fixNegativeConstraintsInUncoarsening(PartitionedHypergraph& partitioned_hg,
                                          const Context& context) {
  LOG << "fixing constraints while uncoarsening";
  gain_cache_t gain_cache = GainCachePtr::constructGainCache(context);
  std::unique_ptr<IRebalancer> rebalancer = RebalancerFactory::getInstance().createObject(
      context.refinement.rebalancing.algorithm, partitioned_hg.initialNumNodes(), context, gain_cache);
  frontToBackConstraints(partitioned_hg, context, gain_cache);
  GainCachePtr::deleteGainCache(gain_cache);
}


template<typename PartitionedHypergraph>
void postprocessNegativeConstraints(PartitionedHypergraph& partitioned_hg,
                                    const Context& context) {
  gain_cache_t gain_cache = GainCachePtr::constructGainCache(context);
  std::unique_ptr<IRebalancer> rebalancer = RebalancerFactory::getInstance().createObject(
      context.refinement.rebalancing.algorithm, partitioned_hg.initialNumNodes(), context, gain_cache);

  LOG << "-------------- stats before postprocessing --------------";
  LOG << "km1       ="<< metrics::quality(partitioned_hg, context);
  LOG << "Imbalance ="<<metrics::imbalance(partitioned_hg, context);
  LOG << (constraintsMet(partitioned_hg)? "Constrains are respected before partitioner" : "! Constrains are not respected before partitioner !");

  frontToBackConstraints(partitioned_hg, context, gain_cache);
  //descendingConstraintDegree(partitioned_hg, context, gain_cache);

  Metrics metrics { metrics::quality(partitioned_hg, context), metrics::imbalance(partitioned_hg, context) };
  mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(partitioned_hg);
  rebalancer->initialize(phg);
  rebalancer->refine(phg, {}, metrics, 0.0);
  LOG << "-------------- stats after postprocessing --------------";
  LOG << "km1       ="<< metrics::quality(partitioned_hg, context);
  LOG << "Imbalance ="<<metrics::imbalance(partitioned_hg, context);
  LOG << (constraintsMet(partitioned_hg)? "Constrains were respected from balancer" : "!!! Balancer destroyed constrains !!!");
  LOG << "";
  GainCachePtr::deleteGainCache(gain_cache);
}

} // namespace mt_kahypar::constraints
/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2017 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2017 Robin Andre <robinandre1995@web.de>
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
#pragma once

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

#include "mt-kahypar/partition/evo_partitioner.h"
#include "mt-kahypar/partition/multilevel.h"
#include "mt-kahypar/partition/evolutionary/edge_frequency.h"
#include "mt-kahypar/partition/evolutionary/population.h"
#include "mt-kahypar/partition/mapping/target_graph.h"
#include "mt-kahypar/partition/partitioner.h"


namespace mt_kahypar::combine {
static constexpr bool debug = false;

  vec<PartitionID> combinePartitions(const std::vector<std::vector<PartitionID>> &parent_partitions);
  Context modifyContext(const Context& context, ContextModifierParameters params);

  template<typename TypeTraits>
  std::vector<PartitionID> createRandomPartition(const typename TypeTraits::Hypergraph& hypergraph, const Context& context) {
    typename TypeTraits::Hypergraph hg = hypergraph.copy(parallel_tag_t{});
    Context c(context);
    typename TypeTraits::PartitionedHypergraph partitioned_hypergraph(context.partition.k, hg);

    const size_t base_seed = context.partition.seed;
    // Randomly assign nodes to blocks
    hg.doParallelForAllNodes([&](const HypernodeID& hn) {
        PartitionID block;
        if (context.partition.deterministic) {
            // derive a per-node seed deterministically
            size_t node_seed = base_seed + static_cast<size_t>(hn);
            std::mt19937 prng(static_cast<uint32_t>(node_seed));
            std::uniform_int_distribution<PartitionID> dist(0, context.partition.k - 1);
            block = dist(prng);
        } else {
            block = utils::Randomize::instance().getRandomInt(0, context.partition.k - 1, THREAD_ID);
        }
        partitioned_hypergraph.setOnlyNodePart(hn, block);
    });

    // Iinitialize partition data structures
    partitioned_hypergraph.initializePartition();

    // V-Cycle to improve partitioning
    Multilevel<TypeTraits>::partitionVCycle(hg, partitioned_hypergraph, c, nullptr);

    // extract partition
    std::vector<PartitionID> partition(hg.initialNumNodes(), 0);
    partitioned_hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
        partition[hn] = partitioned_hypergraph.partID(hn);
    });
    return partition;
  }

  template<typename TypeTraits>
    std::vector<PartitionID> createDegreeSortedPartition(const typename TypeTraits::Hypergraph& hypergraph, const Context& context) {
    std::vector<std::pair<HypernodeID, size_t>> degrees;
    typename TypeTraits::Hypergraph hg = hypergraph.copy(parallel_tag_t{});
    Context c(context);
    typename TypeTraits::PartitionedHypergraph partitioned_hypergraph(context.partition.k, hg);

    for ( const HypernodeID& hn : hg.nodes() ) {
      degrees.emplace_back(hn, hg.nodeDegree(hn));
    }
    std::sort(degrees.begin(), degrees.end(),
              [](const std::pair<HypernodeID, size_t>& a,
                 const std::pair<HypernodeID, size_t>& b) {
                    return a.second > b.second;
              });

    std::vector<PartitionID> partition(hg.initialNumNodes(), 0);

    // Split partition based on degree -- equisize blocks
    size_t block_size = hg.initialNumNodes() / context.partition.k;
    for ( size_t i = 0; i < degrees.size(); ++i ) {
      PartitionID block = std::min(static_cast<PartitionID>(i / block_size), static_cast<PartitionID>(context.partition.k - 1));
      partitioned_hypergraph.setOnlyNodePart(degrees[i].first, block);
    }

    // Iinitialize partition data structures
    partitioned_hypergraph.initializePartition();

    // V-Cycle to improve partitioning
    Multilevel<TypeTraits>::partitionVCycle(hg, partitioned_hypergraph, c, nullptr);

    // extract partition
    partitioned_hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
        partition[hn] = partitioned_hypergraph.partID(hn);
    });
    return partition;

  }
  template<typename TypeTraits>
  Individual runEvolutionVCycle(typename TypeTraits::PartitionedHypergraph& partitioned_hypergraph, typename TypeTraits::Hypergraph& hypergraph, vec<PartitionID>& comms, const Context& context, std::unordered_map<PartitionID, int>& comm_to_block, TargetGraph* target_graph) {
    partitioned_hypergraph.initializePartition();
    hypergraph.setCommunityIDs(std::move(comms));
    if (context.partition.mode == Mode::direct) {
      //V-cycle requires a context with initialized part weights
      Context vc_context(context);
      vc_context.setupPartWeights(hypergraph.totalWeight());
      vec<EdgeMetadata> edge_md;
      Multilevel<TypeTraits>::evolutionPartitionVCycle(hypergraph, partitioned_hypergraph, std::move(edge_md), vc_context, comm_to_block, target_graph);
    } else {
      throw InvalidParameterException("Invalid partitioning mode!");
    }

    return Individual(partitioned_hypergraph, context);
  }

  template<typename TypeTraits>
  Individual combineParentsWithVCycle(const std::vector<std::vector<PartitionID>>& parent_partitions,
                               std::vector<PartitionID> best_partition, const Context& context, TargetGraph* target_graph,
                               const typename TypeTraits::Hypergraph &input_hg) {

    std::unordered_map<PartitionID, int> comm_to_block;
    vec<PartitionID> comms = combinePartitions(parent_partitions);


    typename TypeTraits::Hypergraph hypergraph = input_hg.copy(parallel_tag_t{});
    typename TypeTraits::PartitionedHypergraph partitioned_hypergraph(context.partition.k, hypergraph);

    for ( const HypernodeID& hn : hypergraph.nodes() ) {
      partitioned_hypergraph.setOnlyNodePart(hn, best_partition[hn]);
      if ( comm_to_block.find(comms[hn]) == comm_to_block.end() ) {
        comm_to_block[comms[hn]] = best_partition[hn];
      }
    }
    return runEvolutionVCycle<TypeTraits>(partitioned_hypergraph, hypergraph, comms, context, comm_to_block,
                                          target_graph);
  }

  template <typename TypeTraits>
  Individual usingKWaySelection(const typename TypeTraits::Hypergraph& input_hg, TargetGraph* target_graph, Population& population, const Context& context, std::mt19937* rng) {
    std::vector<size_t> parents;
    //Maybe change to actually use Tournament Selection
    const size_t best(population.sampleKParentsReturnBestIndex(parents, context.evolutionary.kway_combine,
                                                               context.partition.deterministic, rng));

    std::vector<PartitionID> best_partition = population.partitionCopySafe(best);

    // aquire lock --- possibly unnecessary
    std::vector<std::vector<PartitionID>> parent_partitions;
    for (auto parent_id : parents) {
      parent_partitions.push_back(population.partitionCopySafe(parent_id)); // FIXED
    }

    return combineParentsWithVCycle<TypeTraits>(parent_partitions, best_partition, context, target_graph, input_hg);
  }

  template <typename TypeTraits>
  Individual usingArtificialSecondParent(const typename TypeTraits::Hypergraph& input_hg, Population& population, TargetGraph* target_graph, ContextModifierParameters params,const Context& context, std::mt19937* rng) {
      std::vector<std::vector<PartitionID>> parent_partitions;
      size_t best(population.randomIndividualSafe(context.partition.deterministic, rng));

      // generate new parent individual with modified context
      Context modified_context = modifyContext(context, params);

      modified_context.setupPartWeights(input_hg.totalWeight());

      std::vector<PartitionID> modified_partition;

      ASSERT(params.use_random_partitions + params.use_degree_sorted_partitions <= 1,
             "Can only use one of random or degree-sorted partitions");
      if (params.use_random_partitions) {
          modified_partition = combine::createRandomPartition<TypeTraits>(input_hg, modified_context);
      } else if (params.use_degree_sorted_partitions) {
          modified_partition = combine::createDegreeSortedPartition<TypeTraits>(input_hg, modified_context);
      } else {
          modified_partition = EvoPartitioner<TypeTraits>::createPartition(input_hg, modified_context, target_graph);
      }

      std::vector<PartitionID> best_partition = population.partitionCopySafe(best);
      parent_partitions.push_back(best_partition);
      parent_partitions.push_back(modified_partition);

    return combineParentsWithVCycle<TypeTraits>(parent_partitions, best_partition, context, target_graph, input_hg);
  }

  template <typename TypeTraits>
  Individual usingMultiEdgeFrequency(const typename TypeTraits::Hypergraph& input_hg, Population& population,const Context& context, TargetGraph* target_graph, std::mt19937* rng) {
    Context sub_context = Context(context);
    typename TypeTraits::Hypergraph hypergraph = input_hg.copy(parallel_tag_t{});
    size_t parent_amount = std::ceil(std::sqrt(context.evolutionary.population_size));

    if (!sub_context.partition.use_individual_part_weights) {
      sub_context.partition.max_part_weights.clear();
    }

    sub_context.coarsening.algorithm = CoarseningAlgorithm::three_phase_coarsener;
    sub_context.coarsening.rating.degree_similarity_policy = DegreeSimilarityPolicy::guided;

    //parent amount maybe sqrt(population.size) rounded up
    std::vector<size_t> parents;
    population.sampleKParentsReturnBestIndex(parents, parent_amount, context.partition.deterministic, rng);
    //compute edge frequencies
    vec<EdgeMetadata> edge_md(hypergraph.initialNumEdges(), 0);
    for (auto parent_id : parents) {
      std::vector<HyperedgeID> cut_edges = population.cutEdgesCopySave(parent_id);
      for (auto edge : cut_edges) {
        edge_md[edge] += 1.0 / parent_amount;
      }
    }
    //run three phase coarsener with computed edge frequencies of parents
    //individual is created as a complete new partition with new initial partitioning and new coarsening
    typename TypeTraits::PartitionedHypergraph partitioned_hypergraph =
    Partitioner<TypeTraits>::partition(hypergraph, std::move(edge_md), sub_context, target_graph);
    return Individual(partitioned_hypergraph, sub_context);
  }
} // namespace mt_kahypar::combine


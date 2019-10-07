/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

#include "kahypar/meta/policy_registry.h"

#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/preprocessing/single_node_hyperedge_remover.h"
#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/preprocessing/community_redistributor.h"
#include "mt-kahypar/partition/initial_partitioning/initial_partitioner.h"

namespace mt_kahypar {
namespace partition {

class Partitioner {
 private:
  static constexpr bool debug = false;

 public:
  Partitioner() :
    _single_node_he_remover() { }

  Partitioner(const Partitioner&) = delete;
  Partitioner& operator= (const Partitioner&) = delete;

  Partitioner(Partitioner&&) = delete;
  Partitioner& operator= (Partitioner&&) = delete;

  ~Partitioner() = default;

  inline void partition(Hypergraph& hypergraph, Context& context);

 private:

  static inline void setupContext(const Hypergraph& hypergraph, Context& context);

  static inline void configurePreprocessing(const Hypergraph& hypergraph, Context& context);

  inline void sanitize(Hypergraph& hypergraph, const Context& context);

  inline void preprocess(Hypergraph& hypergraph, const Context& context);

  inline void redistribution(Hypergraph& hypergraph, const Context& context);

  inline void postprocess(Hypergraph& hypergraph);

  SingleNodeHyperedgeRemover _single_node_he_remover;
};

inline void Partitioner::setupContext(const Hypergraph& hypergraph, Context& context) {
  context.coarsening.contraction_limit =
    context.coarsening.contraction_limit_multiplier * context.partition.k;

  context.coarsening.hypernode_weight_fraction =
    context.coarsening.max_allowed_weight_multiplier
    / context.coarsening.contraction_limit;

  context.coarsening.max_allowed_node_weight = ceil(context.coarsening.hypernode_weight_fraction
                                                    * hypergraph.totalWeight());

  context.setupPartWeights(hypergraph.totalWeight());

  if ( context.coarsening.use_hypernode_degree_threshold ) {
    // TODO(heuer): replace this with a nice statistical detection of power law distribution
    double avg_hypernode_degree = metrics::avgHypernodeDegree(hypergraph);
    double stdev_hn_degree = 0.0;
    for (const auto& hn : hypergraph.nodes()) {
      stdev_hn_degree += (hypergraph.nodeDegree(hn) - avg_hypernode_degree) *
                        (hypergraph.nodeDegree(hn) - avg_hypernode_degree);
    }
    stdev_hn_degree = std::sqrt(stdev_hn_degree / (hypergraph.currentNumNodes() - 1));
    HyperedgeID rank_hypernode_degree = metrics::hypernodeDegreeRank(hypergraph,
      hypergraph.initialNumNodes() - std::ceil(0.00166 * hypergraph.initialNumNodes()));
    if ( avg_hypernode_degree + 5 * stdev_hn_degree < rank_hypernode_degree && rank_hypernode_degree > 250 ) {
      context.coarsening.hypernode_degree_threshold = rank_hypernode_degree;
    }
  }
}

inline void Partitioner::configurePreprocessing(const Hypergraph& hypergraph, Context& context) {
  unused(hypergraph);
  unused(context);
}

inline void Partitioner::sanitize(Hypergraph& hypergraph, const Context& context) {
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  const auto result = _single_node_he_remover.removeSingleNodeHyperedges(hypergraph);
  if (context.partition.verbose_output && result.num_removed_single_node_hes > 0) {
    LOG << "Performing single-node HE removal:";
    LOG << "\033[1m\033[31m" << " # removed hyperedges with |e|=1 = "
        << result.num_removed_single_node_hes
        << "\033[0m";
    LOG << "\033[1m\033[31m" << " ===>" << result.num_unconnected_hns
        << "unconnected HNs could have been removed" << "\033[0m";
    io::printStripe();
  }
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("single_node_hyperedge_removal", "Single Node Hyperedge Removal",
    "preprocessing", mt_kahypar::utils::Timer::Type::PREPROCESSING, 0, std::chrono::duration<double>(end - start).count());
}

inline void Partitioner::preprocess(Hypergraph& hypergraph, const Context& context) {
  io::printTopLevelPreprocessingBanner(context);

  HighResClockTimepoint global_start = std::chrono::high_resolution_clock::now();
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  std::vector<PartitionID> communities;
  io::readPartitionFile(context.partition.graph_community_filename, communities);
  ASSERT(communities.size() == hypergraph.initialNumNodes());
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("read_community_file", "Read Community File",
    "community_detection", mt_kahypar::utils::Timer::Type::PREPROCESSING, 0, std::chrono::duration<double>(end - start).count());

  // Stream community ids into hypergraph
  start = std::chrono::high_resolution_clock::now();
  tbb::parallel_for(tbb::blocked_range<HypernodeID>(0UL, hypergraph.initialNumNodes()),
    [&](const tbb::blocked_range<HypernodeID>& range) {
    for ( HypernodeID hn = range.begin(); hn < range.end(); ++hn ) {
      hypergraph.streamCommunityID(hypergraph.globalNodeID(hn), communities[hn]);
    }
  });
  end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("stream_community_ids", "Stream Community IDs",
    "community_detection", mt_kahypar::utils::Timer::Type::PREPROCESSING, 1, std::chrono::duration<double>(end - start).count());

  // Initialize Communities
  start = std::chrono::high_resolution_clock::now();
  hypergraph.initializeCommunities();
  end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("initialize_communities", "Initialize Communities",
    "community_detection", mt_kahypar::utils::Timer::Type::PREPROCESSING, 2, std::chrono::duration<double>(end - start).count());
  HighResClockTimepoint global_end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("community_detection", "Community Detection",
    "preprocessing", mt_kahypar::utils::Timer::Type::PREPROCESSING, 1, std::chrono::duration<double>(global_end - global_start).count());

  // Redistribute Hypergraph based on communities
  start = std::chrono::high_resolution_clock::now();
  redistribution(hypergraph, context);
  end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("redistribution", "Redistribution",
    "preprocessing", mt_kahypar::utils::Timer::Type::PREPROCESSING, 2, std::chrono::duration<double>(end - start).count());
}

inline void Partitioner::redistribution(Hypergraph& hypergraph, const Context& context) {
  std::unique_ptr<preprocessing::ICommunityAssignment> community_assignment =
    RedistributionFactory::getInstance().createObject(
      context.shared_memory.assignment_strategy, hypergraph, context);

  std::vector<PartitionID> community_node_mapping = community_assignment->computeAssignment();
  if ( context.shared_memory.use_community_redistribution && TBBNumaArena::instance().num_used_numa_nodes() > 1 ) {
    HyperedgeWeight remote_pin_count_before = metrics::remotePinCount(hypergraph);
    hypergraph = preprocessing::CommunityRedistributor::redistribute(hypergraph, context.partition.k, community_node_mapping);
    HyperedgeWeight remote_pin_count_after = metrics::remotePinCount(hypergraph);
    if ( context.partition.verbose_output ) {
      LOG << "Hypergraph Redistribution Results:";
      LOG << " Remote Pin Count Before Redistribution   =" << remote_pin_count_before;
      LOG << " Remote Pin Count After Redistribution    =" << remote_pin_count_after;
      io::printStripe();
    }
  }
  hypergraph.setCommunityNodeMapping(std::move(community_node_mapping));
}

inline void Partitioner::postprocess(Hypergraph& hypergraph) {
  _single_node_he_remover.restoreSingleNodeHyperedges(hypergraph);
}

inline void Partitioner::partition(Hypergraph& hypergraph, Context& context) {
  configurePreprocessing(hypergraph, context);

  setupContext(hypergraph, context);
  io::printInputInformation(context, hypergraph);

  // ################## PREPROCESSING ##################
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  preprocess(hypergraph, context);
  sanitize(hypergraph, context);
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("preprocessing", "Preprocessing",
    "", mt_kahypar::utils::Timer::Type::PREPROCESSING, 1, std::chrono::duration<double>(end - start).count());

  // ################## COARSENING ##################
  io::printCoarseningBanner(context);
  start = std::chrono::high_resolution_clock::now();
  std::unique_ptr<ICoarsener> coarsener =
    CoarsenerFactory::getInstance().createObject(
      context.coarsening.algorithm, hypergraph, context);
  coarsener->coarsen();
  end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("coarsening", "Coarsening",
    "", mt_kahypar::utils::Timer::Type::COARSENING, 2, std::chrono::duration<double>(end - start).count());

  if ( context.partition.verbose_output ) {
    io::printHypergraphInfo(hypergraph, "Coarsened Hypergraph");
  }

  // ################## INITIAL PARTITIONING ##################
  io::printInitialPartitioningBanner(context);
  start = std::chrono::high_resolution_clock::now();
  InitialPartitioner initial_partitioner(hypergraph, context);
  initial_partitioner.initialPartition();
  end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("initial_partitioning", "Initial Partitioning",
    "", mt_kahypar::utils::Timer::Type::INITIAL_PARTITIONING, 3, std::chrono::duration<double>(end - start).count());

  io::printPartitioningResults(hypergraph, context, "Initial Partitioning Results:");

  // ################## LOCAL SEARCH ##################
  io::printLocalSearchBanner(context);
  start = std::chrono::high_resolution_clock::now();
  coarsener->uncoarsen();
  end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("refinement", "Refinement",
    "", mt_kahypar::utils::Timer::Type::REFINEMENT, 4, std::chrono::duration<double>(end - start).count());

  postprocess(hypergraph);

  if ( context.partition.verbose_output ) {
    io::printHypergraphInfo(hypergraph, "Uncoarsened Hypergraph");
    io::printStripe();
    io::printPartitioningResults(hypergraph, context, "Local Search Results:");
  }
}


} // namespace partition
} // namespace mt_kahypar
/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"

#include "kahypar/meta/policy_registry.h"

#include "mt-kahypar/datastructures/graph.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/initial_partitioning/direct_initial_partitioner.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/multilevel.h"
#include "mt-kahypar/partition/preprocessing/community_detection/parallel_louvain.h"
#include "mt-kahypar/partition/preprocessing/community_reassignment/community_redistributor.h"
#include "mt-kahypar/partition/preprocessing/community_reassignment/single_node_hyperedge_remover.h"
#include "mt-kahypar/utils/stats.h"

namespace mt_kahypar {
namespace partition {
class Partitioner {
 private:
  static constexpr bool debug = false;

 public:
  Partitioner() :
    _single_node_he_remover() { }

  Partitioner(const Partitioner&) = delete;
  Partitioner & operator= (const Partitioner &) = delete;

  Partitioner(Partitioner&&) = delete;
  Partitioner & operator= (Partitioner &&) = delete;

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
  if (context.initial_partitioning.mode == InitialPartitioningMode::direct) {
    context.coarsening.contraction_limit =
      context.coarsening.contraction_limit_multiplier * context.partition.k;
  } else {
    context.coarsening.contraction_limit =
      2 * std::max(context.shared_memory.num_threads, (size_t)context.partition.k) *
      context.coarsening.contraction_limit_multiplier;
  }

  context.coarsening.hypernode_weight_fraction =
    context.coarsening.max_allowed_weight_multiplier
    / context.coarsening.contraction_limit;

  context.coarsening.max_allowed_node_weight = ceil(context.coarsening.hypernode_weight_fraction
                                                    * hypergraph.totalWeight());

  context.setupPartWeights(hypergraph.totalWeight());

  if (context.coarsening.use_hypernode_degree_threshold) {
    // TODO(heuer): replace this with a smarter statistical detection of power law distribution
    double avg_hypernode_degree = metrics::avgHypernodeDegree(hypergraph);
    double stdev_hn_degree = 0.0;
    for (const auto& hn : hypergraph.nodes()) {
      stdev_hn_degree += (hypergraph.nodeDegree(hn) - avg_hypernode_degree) *
                         (hypergraph.nodeDegree(hn) - avg_hypernode_degree);
    }
    stdev_hn_degree = std::sqrt(stdev_hn_degree / (hypergraph.initialNumNodes() - 1));
    HyperedgeID rank_hypernode_degree = metrics::hypernodeDegreeRank(hypergraph,
                                                                     hypergraph.initialNumNodes() - std::ceil(0.00166 * hypergraph.initialNumNodes()));
    if (avg_hypernode_degree + 5 * stdev_hn_degree < rank_hypernode_degree && rank_hypernode_degree > 250) {
      context.coarsening.hypernode_degree_threshold = rank_hypernode_degree;
    }
  }

  context.sanityCheck();
}

inline void Partitioner::configurePreprocessing(const Hypergraph& hypergraph, Context& context) {
  const double density = static_cast<double>(hypergraph.initialNumEdges()) /
                         static_cast<double>(hypergraph.initialNumNodes());
  if (context.preprocessing.community_detection.edge_weight_function == LouvainEdgeWeight::hybrid) {
    if (density < 0.75) {
      context.preprocessing.community_detection.edge_weight_function = LouvainEdgeWeight::degree;
    } else {
      context.preprocessing.community_detection.edge_weight_function = LouvainEdgeWeight::uniform;
    }
  }
}

inline void Partitioner::sanitize(Hypergraph& hypergraph, const Context& context) {
  utils::Timer::instance().start_timer("single_node_hyperedge_removal", "Single Node Hyperedge Removal");
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
  utils::Timer::instance().stop_timer("single_node_hyperedge_removal");
}

inline void Partitioner::preprocess(Hypergraph& hypergraph, const Context& context) {
  io::printTopLevelPreprocessingBanner(context);

  utils::Timer::instance().start_timer("community_detection", "Community Detection");
  utils::Timer::instance().start_timer("perform_community_detection", "Perform Community Detection");
  ds::Clustering communities(0);
  if ( !context.preprocessing.use_community_structure_from_file ) {
    ds::AdjListGraph graph = ds::AdjListStarExpansion::contructGraph(hypergraph, context);
    communities = ParallelModularityLouvain::run(graph, context);   // TODO(lars): give switch for PLM/SLM
    ds::AdjListStarExpansion::restrictClusteringToHypernodes(hypergraph, communities);
  } else {
    io::readPartitionFile(context.partition.graph_community_filename, communities);
  }
  utils::Timer::instance().stop_timer("perform_community_detection");

  // Stream community ids into hypergraph
  utils::Timer::instance().start_timer("stream_community_ids", "Stream Community IDs");
  tbb::parallel_for(tbb::blocked_range<HypernodeID>(0UL, hypergraph.initialNumNodes()),
                    [&](const tbb::blocked_range<HypernodeID>& range) {
        for (HypernodeID hn = range.begin(); hn < range.end(); ++hn) {
          hypergraph.setCommunityID(hypergraph.globalNodeID(hn), communities[hn]);
        }
      });
  utils::Timer::instance().stop_timer("stream_community_ids");

  // Initialize Communities
  utils::Timer::instance().start_timer("initialize_communities", "Initialize Communities");
  hypergraph.initializeCommunities();
  utils::Timer::instance().stop_timer("initialize_communities");

  utils::Stats::instance().add_stat("num_communities", hypergraph.numCommunities());
  utils::Timer::instance().stop_timer("community_detection");

  if ( context.partition.verbose_output ) {
    io::printCommunityInformation(hypergraph);
    io::printStripe();
  }

  // Redistribute Hypergraph based on communities
  utils::Timer::instance().start_timer("redistribution", "Redistribution");
  redistribution(hypergraph, context);
  utils::Timer::instance().stop_timer("redistribution");
}

inline void Partitioner::redistribution(Hypergraph& hypergraph, const Context& context) {
  std::unique_ptr<preprocessing::ICommunityAssignment> community_assignment =
    RedistributionFactory::getInstance().createObject(
      context.preprocessing.community_redistribution.assignment_strategy, hypergraph, context);

  for (int node = 0; node < TBBNumaArena::instance().num_used_numa_nodes(); ++node) {
    utils::Stats::instance().add_stat("initial_hns_on_numa_node_" + std::to_string(node),
                                      (int64_t)hypergraph.initialNumNodes(node));
    utils::Stats::instance().add_stat("initial_hes_on_numa_node_" + std::to_string(node),
                                      (int64_t)hypergraph.initialNumEdges(node));
    utils::Stats::instance().add_stat("initial_pins_on_numa_node_" + std::to_string(node),
                                      (int64_t)hypergraph.initialNumPins(node));
  }

  std::vector<PartitionID> community_node_mapping = community_assignment->computeAssignment();
  if (context.preprocessing.community_redistribution.use_community_redistribution &&
      TBBNumaArena::instance().num_used_numa_nodes() > 1) {
    HyperedgeWeight remote_pin_count_before = metrics::remotePinCount(hypergraph);
    hypergraph = preprocessing::CommunityRedistributor::redistribute(hypergraph, context.partition.k, community_node_mapping);
    HyperedgeWeight remote_pin_count_after = metrics::remotePinCount(hypergraph);
    utils::Stats::instance().add_stat("remote_pin_count_before", remote_pin_count_before);
    utils::Stats::instance().add_stat("remote_pin_count_after", remote_pin_count_after);
    for (int node = 0; node < TBBNumaArena::instance().num_used_numa_nodes(); ++node) {
      utils::Stats::instance().add_stat("hns_on_numa_node_" + std::to_string(node) + "_after_redistribution",
                                        (int64_t)hypergraph.initialNumNodes(node));
      utils::Stats::instance().add_stat("hes_on_numa_node_" + std::to_string(node) + "_after_redistribution",
                                        (int64_t)hypergraph.initialNumEdges(node));
      utils::Stats::instance().add_stat("pins_on_numa_node_" + std::to_string(node) + "_after_redistribution",
                                        (int64_t)hypergraph.initialNumPins(node));
    }
    if (context.partition.verbose_output) {
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
  utils::Timer::instance().start_timer("preprocessing", "Preprocessing");
  preprocess(hypergraph, context);
  sanitize(hypergraph, context);
  utils::Timer::instance().stop_timer("preprocessing");

  // ################## MULTILEVEL ##################
  multilevel::partition(hypergraph, context, true);

  postprocess(hypergraph);

  if (context.partition.verbose_output) {
    io::printHypergraphInfo(hypergraph, "Uncoarsened Hypergraph");
    io::printStripe();
  }
}
}  // namespace partition
}  // namespace mt_kahypar

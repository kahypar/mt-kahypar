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
#include "mt-kahypar/partition/preprocessing/hypergraph_sparsifier.h"
#include "mt-kahypar/partition/preprocessing/community_detection/parallel_louvain.h"
#include "mt-kahypar/partition/preprocessing/community_reassignment/community_redistributor.h"
#include "mt-kahypar/utils/stats.h"

namespace mt_kahypar {
namespace partition {
class Partitioner {
 private:
  static constexpr bool debug = false;

  static double MIN_DEGREE_RANK;
  static double MAX_DEGREE_RANK;
  static HypernodeID STDEV_MAX_DEGREE_THRESHOLD_FACTOR;
  static HypernodeID HIGH_DEGREE_VERTEX_THRESHOLD;

 public:
  Partitioner() :
    _hypergraph_sparsifier() { }

  Partitioner(const Partitioner&) = delete;
  Partitioner & operator= (const Partitioner &) = delete;

  Partitioner(Partitioner&&) = delete;
  Partitioner & operator= (Partitioner &&) = delete;

  inline void partition(Hypergraph& hypergraph, Context& context);

 private:
  static inline void setupContext(Hypergraph& hypergraph, Context& context);

  static inline void configurePreprocessing(const Hypergraph& hypergraph, Context& context);

  inline void sanitize(Hypergraph& hypergraph, const Context& context);

  inline void preprocess(Hypergraph& hypergraph, const Context& context);

  inline void redistribution(Hypergraph& hypergraph, const Context& context);

  inline void postprocess(Hypergraph& hypergraph, const Context& context);

  HypergraphSparsifier _hypergraph_sparsifier;
};

double Partitioner::MIN_DEGREE_RANK = 0.0001;
double Partitioner::MAX_DEGREE_RANK = 0.01;
HypernodeID Partitioner::STDEV_MAX_DEGREE_THRESHOLD_FACTOR = 25;
HypernodeID Partitioner::HIGH_DEGREE_VERTEX_THRESHOLD = 100;

inline void Partitioner::setupContext(Hypergraph& hypergraph, Context& context) {

  context.setupPartWeights(hypergraph.totalWeight());
  context.setupContractionLimit(hypergraph.totalWeight());

  if (context.coarsening.use_high_degree_vertex_threshold) {
    std::vector<HypernodeID> hn_degrees;
    for ( const HypernodeID& hn : hypergraph.nodes() ) {
      hn_degrees.push_back(hypergraph.nodeDegree(hn));
    }
    std::sort(hn_degrees.begin(), hn_degrees.end(),
      [&](const HypernodeID& lhs, const HypernodeID& rhs) {
        return lhs > rhs;
      });
    double last_stdev = 0.0;
    double stdev_at_threshold = 0.0;
    double stdev_threshold = hn_degrees[0] / STDEV_MAX_DEGREE_THRESHOLD_FACTOR;
    HypernodeID prefix_sum = 0;
    HypernodeID prefix_square_sum = 0;
    for ( size_t i = 0; i < MAX_DEGREE_RANK * hn_degrees.size(); ++i ) {
      prefix_sum += hn_degrees[i];
      prefix_square_sum += (hn_degrees[i] * hn_degrees[i]);
      double prefix_avg = static_cast<double>(prefix_sum) / (i + 1);
      // VAR(X) = E(X^2) - avg^2
      double stdev_degree = std::sqrt(static_cast<double>(prefix_square_sum) / (i + 1) - prefix_avg * prefix_avg);

      // We accept the current index i as high degree vertex threshold, if
      //  1.) It is a local maximum of the function stdev(hn_degrees[0:i])
      //  2.) stdev(hn_degrees[0:i]) is greater than a threshold
      //  3.) i is greater than some threshold
      if ( last_stdev > stdev_degree &&
           stdev_degree > stdev_threshold &&
           i > MIN_DEGREE_RANK * hypergraph.initialNumNodes() ) {
        context.coarsening.high_degree_vertex_threshold =
          std::max(hn_degrees[i], HIGH_DEGREE_VERTEX_THRESHOLD);
        break;
      }
      if ( i == MIN_DEGREE_RANK * hypergraph.initialNumNodes() ) {
        stdev_at_threshold = stdev_degree;
      }
      last_stdev = stdev_degree;
    }

    // Fallback, if no high degree vertex threshold was detected,
    // but stdev is still above the threshold
    if ( stdev_at_threshold > stdev_threshold &&
         context.coarsening.high_degree_vertex_threshold == std::numeric_limits<HypernodeID>::max() ) {
        context.coarsening.high_degree_vertex_threshold =
          std::max(hn_degrees[MIN_DEGREE_RANK * hypergraph.initialNumNodes()],
            HIGH_DEGREE_VERTEX_THRESHOLD);
    }

    if ( context.coarsening.high_degree_vertex_threshold != std::numeric_limits<HypernodeID>::max() ) {
      hypergraph.markAllHighDegreeVertices(context.coarsening.high_degree_vertex_threshold);
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
  const HyperedgeID num_removed_single_node_hes =
    _hypergraph_sparsifier.removeSingleNodeHyperedges(hypergraph);
  utils::Timer::instance().stop_timer("single_node_hyperedge_removal");

  utils::Timer::instance().start_timer("degree_zero_hypernode_removal", "Degree Zero Hypernode Removal");
  const HypernodeID num_removed_degree_zero_hypernodes =
    _hypergraph_sparsifier.removeDegreeZeroHypernodes(hypergraph);
  utils::Timer::instance().stop_timer("degree_zero_hypernode_removal");

  if (context.partition.verbose_output && num_removed_single_node_hes > 0) {
    LOG << "Performing single-node HE and degree-zero HN removal:";
    LOG << "\033[1m\033[31m" << " # removed"
        << num_removed_single_node_hes << "hyperedges with |e|=1"
        << "\033[0m";
    LOG << "\033[1m\033[31m" << " # removed"
        << num_removed_degree_zero_hypernodes << "hypernodes with d(v)=0"
        << "\033[0m";
    io::printStripe();
  }
}

inline void Partitioner::preprocess(Hypergraph& hypergraph, const Context& context) {
  io::printTopLevelPreprocessingBanner(context);

  utils::Timer::instance().start_timer("community_detection", "Community Detection");
  utils::Timer::instance().start_timer("perform_community_detection", "Perform Community Detection");
  ds::Clustering communities(0);
  if (!context.preprocessing.use_community_structure_from_file) {
    ds::AdjListGraph graph = ds::AdjListStarExpansion::constructGraph(hypergraph, context, true);
    communities = ParallelModularityLouvain::run(graph, context);   // TODO(lars): give switch for PLM/SLM
    ds::AdjListStarExpansion::restrictClusteringToHypernodes(hypergraph, communities);
    _hypergraph_sparsifier.assignAllDegreeZeroHypernodesToSameCommunity(hypergraph, communities);
  } else {
    io::readPartitionFile(context.partition.graph_community_filename, communities);
  }
  utils::Timer::instance().stop_timer("perform_community_detection");

  // Stream community ids into hypergraph
  utils::Timer::instance().start_timer("stream_community_ids", "Stream Community IDs");
  tbb::parallel_for(tbb::blocked_range<HypernodeID>(0UL, hypergraph.initialNumNodes()),
                    [&](const tbb::blocked_range<HypernodeID>& range) {
        for (HypernodeID id = range.begin(); id < range.end(); ++id) {
          const HypernodeID hn = hypergraph.globalNodeID(id);
          hypergraph.setCommunityID(hn, communities[id]);
        }
      });
  utils::Timer::instance().stop_timer("stream_community_ids");

  // Initialize Communities
  utils::Timer::instance().start_timer("initialize_communities", "Initialize Communities");
  hypergraph.initializeCommunities();
  utils::Timer::instance().stop_timer("initialize_communities");

  utils::Stats::instance().add_stat("num_communities", hypergraph.numCommunities());
  utils::Timer::instance().stop_timer("community_detection");

  if (context.partition.verbose_output) {
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

inline void Partitioner::postprocess(Hypergraph& hypergraph, const Context& context) {
  _hypergraph_sparsifier.restoreDegreeZeroHypernodes(hypergraph, context);
  _hypergraph_sparsifier.restoreSingleNodeHyperedges(hypergraph);
}

inline void Partitioner::partition(Hypergraph& hypergraph, Context& context) {
  configurePreprocessing(hypergraph, context);
  setupContext(hypergraph, context);

  io::printContext(context);
  io::printInputInformation(context, hypergraph);

  // ################## PREPROCESSING ##################
  utils::Timer::instance().start_timer("preprocessing", "Preprocessing");
  preprocess(hypergraph, context);
  sanitize(hypergraph, context);
  utils::Timer::instance().stop_timer("preprocessing");

  // ################## MULTILEVEL ##################
  multilevel::partition(hypergraph, context, true, TBBNumaArena::instance());

  postprocess(hypergraph, context);

  if (context.partition.verbose_output) {
    io::printHypergraphInfo(hypergraph, "Uncoarsened Hypergraph");
    io::printStripe();
  }
}
}  // namespace partition
}  // namespace mt_kahypar

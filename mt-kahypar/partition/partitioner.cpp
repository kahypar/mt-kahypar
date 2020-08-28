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

#include "partitioner.h"

#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/partition/multilevel.h"
#include "mt-kahypar/partition/preprocessing/sparsification/degree_zero_hn_remover.h"
#include "mt-kahypar/partition/preprocessing/sparsification/large_he_remover.h"
#include "mt-kahypar/partition/preprocessing/community_detection/parallel_louvain.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"



namespace mt_kahypar {

  void setupContext(Hypergraph& hypergraph, Context& context) {
    context.partition.large_hyperedge_size_threshold = std::max(hypergraph.initialNumNodes() *
                                                                context.partition.large_hyperedge_size_threshold_factor, 100.0);
    context.setupPartWeights(hypergraph.totalWeight());
    context.setupContractionLimit(hypergraph.totalWeight());
    context.sanityCheck();
  }

  void configurePreprocessing(const Hypergraph& hypergraph, Context& context) {
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

  void sanitize(Hypergraph& hypergraph, Context& context,
                DegreeZeroHypernodeRemover& degree_zero_hn_remover,
                LargeHyperedgeRemover& large_he_remover) {

    utils::Timer::instance().start_timer("degree_zero_hypernode_removal", "Degree Zero Hypernode Removal");
    const HypernodeID num_removed_degree_zero_hypernodes =
            degree_zero_hn_remover.removeDegreeZeroHypernodes(hypergraph);
    utils::Timer::instance().stop_timer("degree_zero_hypernode_removal");

    utils::Timer::instance().start_timer("large_hyperedge_removal", "Large Hyperedge Removal");
    const HypernodeID num_removed_large_hyperedges =
            large_he_remover.removeLargeHyperedges(hypergraph);
    utils::Timer::instance().stop_timer("large_hyperedge_removal");

    const HyperedgeID num_removed_single_node_hes = hypergraph.numRemovedHyperedges();
    if (context.partition.verbose_output &&
        ( num_removed_single_node_hes > 0 ||
          num_removed_degree_zero_hypernodes > 0 ||
          num_removed_large_hyperedges > 0 )) {
      LOG << "Performed single-node/large HE removal and degree-zero HN contractions:";
      LOG << "\033[1m\033[31m" << " # removed"
          << num_removed_single_node_hes << "single-pin hyperedges during hypergraph file parsing"
          << "\033[0m";
      LOG << "\033[1m\033[31m" << " # removed"
          << num_removed_large_hyperedges << "large hyperedges with |e| >" << large_he_remover.largeHyperedgeThreshold() << "\033[0m";
      LOG << "\033[1m\033[31m" << " # contracted"
          << num_removed_degree_zero_hypernodes << "hypernodes with d(v) = 0 and w(v) = 1"
          << "\033[0m";
      io::printStripe();
    }
  }

  void preprocess(Hypergraph& hypergraph, Context& context) {
    if ( context.preprocessing.use_community_detection ) {
      io::printTopLevelPreprocessingBanner(context);

      utils::Timer::instance().start_timer("community_detection", "Community Detection");
      utils::Timer::instance().start_timer("construct_graph", "Construct Graph");
      Graph graph(hypergraph, context.preprocessing.community_detection.edge_weight_function);
      utils::Timer::instance().stop_timer("construct_graph");
      utils::Timer::instance().start_timer("perform_community_detection", "Perform Community Detection");
      ds::Clustering communities = community_detection::run_parallel_louvain(graph, context);
      graph.restrictClusteringToHypernodes(hypergraph, communities);
      hypergraph.setCommunityIDs(std::move(communities));
      utils::Timer::instance().stop_timer("perform_community_detection");
      utils::Timer::instance().stop_timer("community_detection");

      if (context.partition.verbose_output) {
        io::printCommunityInformation(hypergraph);
      }
    }
    parallel::MemoryPool::instance().release_mem_group("Preprocessing");
  }


  PartitionedHypergraph partitionVCycle(Hypergraph& hypergraph,
                                        PartitionedHypergraph&& partitioned_hypergraph,
                                        Context& context) {
    ASSERT(context.partition.num_vcycles > 0);

    for ( size_t i = 0; i < context.partition.num_vcycles; ++i ) {
      // Reset memory pool
      hypergraph.reset();
      parallel::MemoryPool::instance().reset();
      parallel::MemoryPool::instance().release_mem_group("Preprocessing");

      // Store partition and assign it as community ids in order to
      // restrict contractions in v-cycle to partition ids
      hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
        hypergraph.setCommunityID(hn, partitioned_hypergraph.partID(hn));
      });

      // V-Cycle Multilevel Partitioning
      io::printVCycleBanner(context, i + 1);
      // TODO why does this have to make a copy
      partitioned_hypergraph = multilevel::partition(
              hypergraph, context, true, TBBNumaArena::GLOBAL_TASK_GROUP, true /* vcycle */);
    }

    return std::move(partitioned_hypergraph);
  }

  PartitionedHypergraph partition(Hypergraph& hypergraph, Context& context) {
    configurePreprocessing(hypergraph, context);
    setupContext(hypergraph, context);

    io::printContext(context);
    io::printMemoryPoolConsumption(context);
    io::printInputInformation(context, hypergraph);

    // ################## PREPROCESSING ##################
    utils::Timer::instance().start_timer("preprocessing", "Preprocessing");
    preprocess(hypergraph, context);

    DegreeZeroHypernodeRemover degree_zero_hn_remover(context);
    LargeHyperedgeRemover large_he_remover(context);
    sanitize(hypergraph, context, degree_zero_hn_remover, large_he_remover);
    utils::Timer::instance().stop_timer("preprocessing");

    // ################## MULTILEVEL ##################
    PartitionedHypergraph partitioned_hypergraph = multilevel::partition(
            hypergraph, context, true, TBBNumaArena::GLOBAL_TASK_GROUP);

    // ################## V-Cycle s##################
    if ( context.partition.num_vcycles > 0 ) {
      partitioned_hypergraph = partitionVCycle(hypergraph, std::move(partitioned_hypergraph), context);
    }

    // ################## POSTPROCESSING ##################
    utils::Timer::instance().start_timer("postprocessing", "Postprocessing");
    large_he_remover.restoreLargeHyperedges(partitioned_hypergraph);
    degree_zero_hn_remover.restoreDegreeZeroHypernodes(partitioned_hypergraph);
    utils::Timer::instance().stop_timer("postprocessing");

    if (context.partition.verbose_output) {
      io::printHypergraphInfo(partitioned_hypergraph.hypergraph(), "Uncoarsened Hypergraph",
                              context.partition.show_memory_consumption);
      io::printStripe();
    }

    return partitioned_hypergraph;
  }



}
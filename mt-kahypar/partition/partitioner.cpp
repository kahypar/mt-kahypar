/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "partitioner.h"

#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/partition/multilevel.h"
#include "mt-kahypar/partition/preprocessing/sparsification/degree_zero_hn_remover.h"
#include "mt-kahypar/partition/preprocessing/sparsification/large_he_remover.h"
#include "mt-kahypar/partition/preprocessing/community_detection/parallel_louvain.h"
#include "mt-kahypar/partition/recursive_bipartitioning.h"
#include "mt-kahypar/partition/deep_multilevel.h"
#include "mt-kahypar/utils/hypergraph_statistics.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"


namespace mt_kahypar {

  void setupContext(Hypergraph& hypergraph, Context& context) {
    context.partition.large_hyperedge_size_threshold = std::max(hypergraph.initialNumNodes() *
                                                                context.partition.large_hyperedge_size_threshold_factor, 100.0);
    context.sanityCheck();
    context.setupPartWeights(hypergraph.totalWeight());
    context.setupContractionLimit(hypergraph.totalWeight());
    context.setupThreadsPerFlowSearch();

    // Setup enabled IP algorithms
    if ( context.initial_partitioning.enabled_ip_algos.size() > 0 &&
         context.initial_partitioning.enabled_ip_algos.size() <
         static_cast<size_t>(InitialPartitioningAlgorithm::UNDEFINED) ) {
      ERROR("Size of enabled IP algorithms vector is smaller than number of IP algorithms!");
    } else if ( context.initial_partitioning.enabled_ip_algos.size() == 0 ) {
      context.initial_partitioning.enabled_ip_algos.assign(
        static_cast<size_t>(InitialPartitioningAlgorithm::UNDEFINED), true);
    } else {
      bool is_one_ip_algo_enabled = false;
      for ( size_t i = 0; i < context.initial_partitioning.enabled_ip_algos.size(); ++i ) {
        is_one_ip_algo_enabled |= context.initial_partitioning.enabled_ip_algos[i];
      }
      if ( !is_one_ip_algo_enabled ) {
        ERROR("At least one initial partitioning algorithm must be enabled!");
      }
    }
  }

  void configurePreprocessing(const Hypergraph& hypergraph, Context& context) {
    const double density = static_cast<double>(Hypergraph::is_graph ? hypergraph.initialNumEdges() / 2 : hypergraph.initialNumEdges()) /
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

    utils::Timer& timer = utils::Utilities::instance().getTimer(context.utility_id);
    timer.start_timer("degree_zero_hypernode_removal", "Degree Zero Hypernode Removal");
    const HypernodeID num_removed_degree_zero_hypernodes =
            degree_zero_hn_remover.removeDegreeZeroHypernodes(hypergraph);
    timer.stop_timer("degree_zero_hypernode_removal");

    timer.start_timer("large_hyperedge_removal", "Large Hyperedge Removal");
    const HypernodeID num_removed_large_hyperedges =
            large_he_remover.removeLargeHyperedges(hypergraph);
    timer.stop_timer("large_hyperedge_removal");

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

  bool isGraph(const Hypergraph& hypergraph) {
    if (Hypergraph::is_graph) {
      return true;
    }
    return tbb::parallel_reduce(tbb::blocked_range<HyperedgeID>(
            ID(0), hypergraph.initialNumEdges()), true, [&](const tbb::blocked_range<HyperedgeID>& range, bool isGraph) {
      if ( isGraph ) {
        bool tmp_is_graph = isGraph;
        for (HyperedgeID he = range.begin(); he < range.end(); ++he) {
          if ( hypergraph.edgeIsEnabled(he) ) {
            tmp_is_graph &= (hypergraph.edgeSize(he) == 2);
          }
        }
        return tmp_is_graph;
      }
      return false;
    }, [&](const bool lhs, const bool rhs) {
      return lhs && rhs;
    });
  }

  bool isMeshGraph(const Hypergraph& graph) {
    const HypernodeID num_nodes = graph.initialNumNodes();
    const double avg_hn_degree = utils::avgHypernodeDegree(graph);
    std::vector<HyperedgeID> hn_degrees;
    hn_degrees.resize(graph.initialNumNodes());
    graph.doParallelForAllNodes([&](const HypernodeID& hn) {
      hn_degrees[hn] = graph.nodeDegree(hn);
    });
    const double stdev_hn_degree = utils::parallel_stdev(hn_degrees, avg_hn_degree, num_nodes);
    if (stdev_hn_degree > avg_hn_degree / 2) {
      return false;
    }

    // test whether 99.9th percentile hypernode degree is at most 4 times the average degree
    tbb::enumerable_thread_specific<size_t> num_high_degree_nodes(0);
    graph.doParallelForAllNodes([&](const HypernodeID& node) {
      if (graph.nodeDegree(node) > 4 * avg_hn_degree) {
        num_high_degree_nodes.local() += 1;
      }
    });
    return num_high_degree_nodes.combine(std::plus<>()) <= num_nodes / 1000;
  }

  void preprocess(Hypergraph& hypergraph, Context& context) {
    bool use_community_detection = context.preprocessing.use_community_detection;
    bool is_graph = false;

    utils::Timer& timer = utils::Utilities::instance().getTimer(context.utility_id);
    if ( context.preprocessing.use_community_detection ) {
      timer.start_timer("detect_graph_structure", "Detect Graph Structure");
      is_graph = isGraph(hypergraph);
      if ( is_graph && context.preprocessing.disable_community_detection_for_mesh_graphs ) {
        use_community_detection = !isMeshGraph(hypergraph);
      }
      timer.stop_timer("detect_graph_structure");
    }

    if ( use_community_detection ) {
      io::printTopLevelPreprocessingBanner(context);

      timer.start_timer("community_detection", "Community Detection");
      timer.start_timer("construct_graph", "Construct Graph");
      Graph graph(hypergraph, context.preprocessing.community_detection.edge_weight_function, is_graph);
      if ( !context.preprocessing.community_detection.low_memory_contraction ) {
        graph.allocateContractionBuffers();
      }
      timer.stop_timer("construct_graph");
      timer.start_timer("perform_community_detection", "Perform Community Detection");
      ds::Clustering communities = community_detection::run_parallel_louvain(graph, context);
      graph.restrictClusteringToHypernodes(hypergraph, communities);
      hypergraph.setCommunityIDs(std::move(communities));
      timer.stop_timer("perform_community_detection");
      timer.stop_timer("community_detection");

      if (context.partition.verbose_output) {
        io::printCommunityInformation(hypergraph);
      }
    }
    parallel::MemoryPool::instance().release_mem_group("Preprocessing");
  }

  PartitionedHypergraph partition(Hypergraph& hypergraph, Context& context) {
    configurePreprocessing(hypergraph, context);
    setupContext(hypergraph, context);

    io::printContext(context);
    io::printMemoryPoolConsumption(context);
    io::printInputInformation(context, hypergraph);

    // ################## PREPROCESSING ##################
    utils::Timer& timer = utils::Utilities::instance().getTimer(context.utility_id);
    timer.start_timer("preprocessing", "Preprocessing");
    preprocess(hypergraph, context);

    DegreeZeroHypernodeRemover degree_zero_hn_remover(context);
    LargeHyperedgeRemover large_he_remover(context);
    sanitize(hypergraph, context, degree_zero_hn_remover, large_he_remover);
    timer.stop_timer("preprocessing");

    // ################## MULTILEVEL & VCYCLE ##################
    PartitionedHypergraph partitioned_hypergraph;
    if (context.partition.mode == Mode::direct) {
      partitioned_hypergraph = multilevel::partition(hypergraph, context);
    } else if (context.partition.mode == Mode::recursive_bipartitioning) {
      partitioned_hypergraph = recursive_bipartitioning::partition(hypergraph, context);
    } else if (context.partition.mode == Mode::deep_multilevel) {
      partitioned_hypergraph = deep_multilevel::partition(hypergraph, context);
    } else {
      ERROR("Invalid mode: " << context.partition.mode);
    }

    // ################## POSTPROCESSING ##################
    timer.start_timer("postprocessing", "Postprocessing");
    large_he_remover.restoreLargeHyperedges(partitioned_hypergraph);
    degree_zero_hn_remover.restoreDegreeZeroHypernodes(partitioned_hypergraph);
    timer.stop_timer("postprocessing");

    if (context.partition.verbose_output) {
      io::printHypergraphInfo(partitioned_hypergraph.hypergraph(), "Uncoarsened Hypergraph",
                              context.partition.show_memory_consumption);
      io::printStripe();
    }

    return partitioned_hypergraph;
  }


  void partitionVCycle(PartitionedHypergraph& partitioned_hg, Context& context) {
    Hypergraph& hypergraph = partitioned_hg.hypergraph();
    configurePreprocessing(hypergraph, context);
    setupContext(hypergraph, context);

    io::printContext(context);
    io::printMemoryPoolConsumption(context);
    io::printInputInformation(context, hypergraph);
    io::printPartitioningResults(partitioned_hg, context, "\nInput Partition:");

    // ################## PREPROCESSING ##################
    utils::Timer& timer = utils::Utilities::instance().getTimer(context.utility_id);
    timer.start_timer("preprocessing", "Preprocessing");
    DegreeZeroHypernodeRemover degree_zero_hn_remover(context);
    LargeHyperedgeRemover large_he_remover(context);
    sanitize(hypergraph, context, degree_zero_hn_remover, large_he_remover);
    timer.stop_timer("preprocessing");

    // ################## MULTILEVEL & VCYCLE ##################
    if (context.partition.mode == Mode::direct) {
      multilevel::partitionVCycle(hypergraph, partitioned_hg, context);
    } else {
      ERROR("Invalid V-cycle mode: " << context.partition.mode);
    }

    // ################## POSTPROCESSING ##################
    timer.start_timer("postprocessing", "Postprocessing");
    large_he_remover.restoreLargeHyperedges(partitioned_hg);
    degree_zero_hn_remover.restoreDegreeZeroHypernodes(partitioned_hg);
    timer.stop_timer("postprocessing");

    if (context.partition.verbose_output) {
      io::printHypergraphInfo(partitioned_hg.hypergraph(),
        "Uncoarsened Hypergraph", context.partition.show_memory_consumption);
      io::printStripe();
    }
  }



}
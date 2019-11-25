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

#include <libkahypar.h>

#include "kahypar/partition/context.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/metrics.h"

namespace mt_kahypar {
struct KaHyParPartitioningResult {
  KaHyParPartitioningResult() :
    objective(0),
    imbalance(1.0),
    partition() { }

  explicit KaHyParPartitioningResult(const kahypar_hypernode_id_t num_vertices) :
    objective(0),
    imbalance(1.0),
    partition(num_vertices, -1) { }

  kahypar_hyperedge_weight_t objective;
  double imbalance;
  parallel::scalable_vector<kahypar_partition_id_t> partition;
};

struct KaHyParHypergraph {
  KaHyParHypergraph() :
    num_vertices(0),
    num_hyperedges(0),
    hyperedge_indices(),
    hyperedges(),
    vertex_weights(),
    hyperedge_weights(),
    reverse_mapping(),
    part_weights() { }

  void computePartWeights(const PartitionID k,
                          const parallel::scalable_vector<kahypar_partition_id_t>& partition) {
    ASSERT(partition.size() == vertex_weights.size());
    part_weights.assign(k, 0);
    for (HypernodeID hn = 0; hn < partition.size(); ++hn) {
      PartitionID part = partition[hn];
      ASSERT(part != -1 && part < k);
      part_weights[part] += vertex_weights[hn];
    }
  }

  HypernodeWeight partWeight(const PartitionID id) const {
    ASSERT(id < (PartitionID)part_weights.size());
    return part_weights[id];
  }

  kahypar_hypernode_id_t num_vertices;
  kahypar_hyperedge_id_t num_hyperedges;
  parallel::scalable_vector<size_t> hyperedge_indices;
  parallel::scalable_vector<kahypar_hyperedge_id_t> hyperedges;
  parallel::scalable_vector<kahypar_hypernode_weight_t> vertex_weights;
  parallel::scalable_vector<kahypar_hyperedge_weight_t> hyperedge_weights;
  parallel::scalable_vector<HypernodeID> reverse_mapping;
  parallel::scalable_vector<HypernodeWeight> part_weights;
};

// ! Converts the MT-KaHyPar Hypergraph into raw hypergraph data format taken
// ! by the KaHyPar-Library Interface.
template <typename HyperGraph>
static KaHyParHypergraph convertToKaHyParHypergraph(const HyperGraph& hypergraph,
                                                    const std::vector<HypernodeID>& node_mapping) {
  ASSERT(hypergraph.initialNumNodes() == node_mapping.size());
  KaHyParHypergraph kahypar_hypergraph;

  // Setup hypernodes
  for (const HypernodeID& hn : hypergraph.nodes()) {
    kahypar_hypergraph.reverse_mapping.emplace_back(hn);
    kahypar_hypergraph.vertex_weights.emplace_back(hypergraph.nodeWeight(hn));
    ++kahypar_hypergraph.num_vertices;
  }

  // Setup hyperedges
  kahypar_hypergraph.hyperedge_indices.emplace_back(0);
  for (const HyperedgeID& he : hypergraph.edges()) {
    ++kahypar_hypergraph.num_hyperedges;
    for (const HypernodeID& hn : hypergraph.pins(he)) {
      ASSERT(hypergraph.originalNodeID(hn) < hypergraph.initialNumNodes());
      ASSERT(node_mapping[hypergraph.originalNodeID(hn)] != std::numeric_limits<HypernodeID>::max());
      kahypar_hypergraph.hyperedges.emplace_back(node_mapping[hypergraph.originalNodeID(hn)]);
    }
    kahypar_hypergraph.hyperedge_indices.emplace_back(kahypar_hypergraph.hyperedges.size());
    kahypar_hypergraph.hyperedge_weights.emplace_back(hypergraph.edgeWeight(he));
  }

  return kahypar_hypergraph;
}

// ! Extracts a block of a MT-KaHyPar Hypergraph as a raw hypergraph data format taken
// ! by the KaHyPar-Library Interface.
template <typename HyperGraph>
static KaHyParHypergraph extractBlockAsKaHyParHypergraph(const HyperGraph& hypergraph,
                                                         const PartitionID block,
                                                         std::vector<HypernodeID>& node_mapping,
                                                         const bool cut_net_splitting = true) {
  ASSERT(hypergraph.initialNumNodes() == node_mapping.size());
  ASSERT(block >= 0 && block < hypergraph.k());
  KaHyParHypergraph kahypar_hypergraph;

  // Setup hypernodes
  for (const HypernodeID& hn : hypergraph.nodes()) {
    ASSERT(hypergraph.partID(hn) != -1);
    if (hypergraph.partID(hn) == block) {
      node_mapping[hypergraph.originalNodeID(hn)] = kahypar_hypergraph.num_vertices++;
      kahypar_hypergraph.reverse_mapping.emplace_back(hn);
      kahypar_hypergraph.vertex_weights.emplace_back(hypergraph.nodeWeight(hn));
    }
  }

  // Setup hyperedges
  kahypar_hypergraph.hyperedge_indices.emplace_back(0);
  for (const HyperedgeID& he : hypergraph.edges()) {
    if (hypergraph.pinCountInPart(he, block) > 0 &&
        (hypergraph.connectivity(he) == 1 || cut_net_splitting)) {
      ++kahypar_hypergraph.num_hyperedges;
      for (const HypernodeID& hn : hypergraph.pins(he)) {
        if (hypergraph.partID(hn) == block) {
          ASSERT(hypergraph.originalNodeID(hn) < hypergraph.initialNumNodes());
          ASSERT(node_mapping[hypergraph.originalNodeID(hn)] != std::numeric_limits<HypernodeID>::max());
          ASSERT(hypergraph.partID(hn) != -1);
          kahypar_hypergraph.hyperedges.emplace_back(node_mapping[hypergraph.originalNodeID(hn)]);
        }
      }
      kahypar_hypergraph.hyperedge_indices.emplace_back(kahypar_hypergraph.hyperedges.size());
      kahypar_hypergraph.hyperedge_weights.emplace_back(hypergraph.edgeWeight(he));
    }
  }

  return kahypar_hypergraph;
}

static void kahypar_partition(const KaHyParHypergraph& hypergraph,
                              kahypar::Context& context,
                              KaHyParPartitioningResult& partition) {
  kahypar_partition(hypergraph.num_vertices,
                    hypergraph.num_hyperedges,
                    context.partition.epsilon,
                    context.partition.k,
                    hypergraph.vertex_weights.data(),
                    hypergraph.hyperedge_weights.data(),
                    hypergraph.hyperedge_indices.data(),
                    hypergraph.hyperedges.data(),
                    &partition.objective,
                    reinterpret_cast<kahypar_context_t*>(&context),
                    partition.partition.data());
}

// ! Computes the imbalance of a partition without applying it to the hypergraph.
static double imbalance(KaHyParHypergraph& kahypar_hypergraph,
                        const Context& context,
                        const KaHyParPartitioningResult& kahypar_result) {
  kahypar_hypergraph.computePartWeights(context.partition.k, kahypar_result.partition);
  return metrics::imbalance(kahypar_hypergraph, context);
}

static void sanitizeCheck(kahypar::Context& context) {
  if (context.partition.objective == kahypar::Objective::km1) {
    switch (context.local_search.algorithm) {
      case kahypar::RefinementAlgorithm::kway_fm:
        context.local_search.algorithm = kahypar::RefinementAlgorithm::kway_fm_km1;
        WARNING("Auto switch refinement algorithm from kway_fm to kway_fm_km1 (KaHyPar IP)");
        break;
      case kahypar::RefinementAlgorithm::kway_fm_flow:
        context.local_search.algorithm = kahypar::RefinementAlgorithm::kway_fm_flow_km1;
        WARNING("Auto switch refinement algorithm from kway_fm_flow to kway_fm_flow_km1 (KaHyPar IP)");
        break;
      default:
        break;
    }

    if (context.partition.k > 2) {
      switch (context.local_search.algorithm) {
        case kahypar::RefinementAlgorithm::twoway_flow:
          context.local_search.algorithm = kahypar::RefinementAlgorithm::kway_flow;
          WARNING("Auto switch refinement algorithm from twoway_flow to kway_flow (KaHyPar IP)");
          break;
        case kahypar::RefinementAlgorithm::twoway_fm:
          context.local_search.algorithm = kahypar::RefinementAlgorithm::kway_fm_km1;
          WARNING("Auto switch refinement algorithm from twoway_fm to kway_fm_km1 (KaHyPar IP)");
          break;
        case kahypar::RefinementAlgorithm::twoway_fm_flow:
          context.local_search.algorithm = kahypar::RefinementAlgorithm::kway_fm_flow_km1;
          WARNING("Auto switch refinement algorithm from twoway_fm_flow to kway_fm_flow_km1 (KaHyPar IP)");
          break;
        default:
          break;
      }
    }
  } else {
    switch (context.local_search.algorithm) {
      case kahypar::RefinementAlgorithm::kway_fm_km1:
        context.local_search.algorithm = kahypar::RefinementAlgorithm::kway_fm;
        WARNING("Auto switch refinement algorithm from kway_fm_km1 to kway_fm (KaHyPar IP)");
        break;
      case kahypar::RefinementAlgorithm::kway_fm_flow_km1:
        context.local_search.algorithm = kahypar::RefinementAlgorithm::kway_fm_flow;
        WARNING("Auto switch refinement algorithm from kway_fm_flow_km1 to kway_fm_flow (KaHyPar IP)");
        break;
      default:
        break;
    }

    if (context.partition.k > 2) {
      switch (context.local_search.algorithm) {
        case kahypar::RefinementAlgorithm::twoway_flow:
          context.local_search.algorithm = kahypar::RefinementAlgorithm::kway_flow;
          WARNING("Auto switch refinement algorithm from twoway_flow to kway_flow (KaHyPar IP)");
          break;
        case kahypar::RefinementAlgorithm::twoway_fm:
          context.local_search.algorithm = kahypar::RefinementAlgorithm::kway_fm;
          WARNING("Auto switch refinement algorithm from twoway_fm to kway_fm (KaHyPar IP)");
          break;
        case kahypar::RefinementAlgorithm::twoway_fm_flow:
          context.local_search.algorithm = kahypar::RefinementAlgorithm::kway_fm_flow;
          WARNING("Auto switch refinement algorithm from twoway_fm_flow to kway_fm_flow (KaHyPar IP)");
          break;
        default:
          break;
      }
    }
  }
}

// ! Reads a KaHyPar context file
static kahypar_context_t* readContext(const std::string& context_file) {
  kahypar_context_t* context = kahypar_context_new();
  kahypar_configure_context_from_file(context, context_file.c_str());
  return context;
}

static kahypar_context_t* setupContext(const Context& context, const bool debug) {
  kahypar_context_t* kahypar_c = readContext(context.initial_partitioning.context_file);
  kahypar::Context& kahypar_context = *reinterpret_cast<kahypar::Context*>(kahypar_c);
  kahypar_context.partition.objective = context.partition.objective;
  kahypar_context.partition.epsilon = context.partition.epsilon;
  kahypar_context.partition.k = context.partition.k;
  kahypar_context.partition.seed = context.partition.seed;
  kahypar_context.preprocessing.enable_deduplication = true;
  kahypar_context.preprocessing.enable_min_hash_sparsifier = false;
  kahypar_context.preprocessing.enable_community_detection = false;
  kahypar_context.partition.verbose_output = debug;
  sanitizeCheck(kahypar_context);
  return kahypar_c;
}
}  // namespace mt_kahypar

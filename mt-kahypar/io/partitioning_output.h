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

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/memory_tree.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
namespace io {
namespace internal {
struct Statistic {
  uint64_t min = 0;
  uint64_t q1 = 0;
  uint64_t med = 0;
  uint64_t q3 = 0;
  uint64_t max = 0;
  double avg = 0.0;
  double sd = 0.0;
};

template <typename T>
Statistic createStats(const std::vector<T>& vec, const double avg, const double stdev) {
  internal::Statistic stats;
  if (!vec.empty()) {
    const auto quartiles = kahypar::math::firstAndThirdQuartile(vec);
    stats.min = vec.front();
    stats.q1 = quartiles.first;
    stats.med = kahypar::math::median(vec);
    stats.q3 = quartiles.second;
    stats.max = vec.back();
    stats.avg = avg;
    stats.sd = stdev;
  }
  return stats;
}

void printHypergraphStats(const Statistic& he_size_stats,
                          const Statistic& he_weight_stats,
                          const Statistic& hn_deg_stats,
                          const Statistic& hn_weight_stats) {
  // default double precision is 7
  const uint8_t double_width = 7;
  const uint8_t he_size_width = std::max(kahypar::math::digits(he_size_stats.max), double_width) + 4;
  const uint8_t he_weight_width = std::max(kahypar::math::digits(he_weight_stats.max), double_width) + 4;
  const uint8_t hn_deg_width = std::max(kahypar::math::digits(hn_deg_stats.max), double_width) + 4;
  const uint8_t hn_weight_width = std::max(kahypar::math::digits(hn_weight_stats.max), double_width) + 4;

  LOG << "HE size" << std::right << std::setw(he_size_width + 10)
      << "HE weight" << std::right << std::setw(he_weight_width + 8)
      << "HN degree" << std::right << std::setw(hn_deg_width + 8)
      << "HN weight";
  LOG << "| min=" << std::left << std::setw(he_size_width) << he_size_stats.min
      << " | min=" << std::left << std::setw(he_weight_width) << he_weight_stats.min
      << " | min=" << std::left << std::setw(hn_deg_width) << hn_deg_stats.min
      << " | min=" << std::left << std::setw(hn_weight_width) << hn_weight_stats.min;
  LOG << "| Q1 =" << std::left << std::setw(he_size_width) << he_size_stats.q1
      << " | Q1 =" << std::left << std::setw(he_weight_width) << he_weight_stats.q1
      << " | Q1 =" << std::left << std::setw(hn_deg_width) << hn_deg_stats.q1
      << " | Q1 =" << std::left << std::setw(hn_weight_width) << hn_weight_stats.q1;
  LOG << "| med=" << std::left << std::setw(he_size_width) << he_size_stats.med
      << " | med=" << std::left << std::setw(he_weight_width) << he_weight_stats.med
      << " | med=" << std::left << std::setw(hn_deg_width) << hn_deg_stats.med
      << " | med=" << std::left << std::setw(hn_weight_width) << hn_weight_stats.med;
  LOG << "| Q3 =" << std::left << std::setw(he_size_width) << he_size_stats.q3
      << " | Q3 =" << std::left << std::setw(he_weight_width) << he_weight_stats.q3
      << " | Q3 =" << std::left << std::setw(hn_deg_width) << hn_deg_stats.q3
      << " | Q3 =" << std::left << std::setw(hn_weight_width) << hn_weight_stats.q3;
  LOG << "| max=" << std::left << std::setw(he_size_width) << he_size_stats.max
      << " | max=" << std::left << std::setw(he_weight_width) << he_weight_stats.max
      << " | max=" << std::left << std::setw(hn_deg_width) << hn_deg_stats.max
      << " | max=" << std::left << std::setw(hn_weight_width) << hn_weight_stats.max;
  LOG << "| avg=" << std::left << std::setw(he_size_width) << he_size_stats.avg
      << " | avg=" << std::left << std::setw(he_weight_width) << he_weight_stats.avg
      << " | avg=" << std::left << std::setw(hn_deg_width) << hn_deg_stats.avg
      << " | avg=" << std::left << std::setw(hn_weight_width) << hn_weight_stats.avg;
  LOG << "| sd =" << std::left << std::setw(he_size_width) << he_size_stats.sd
      << " | sd =" << std::left << std::setw(he_weight_width) << he_weight_stats.sd
      << " | sd =" << std::left << std::setw(hn_deg_width) << hn_deg_stats.sd
      << " | sd =" << std::left << std::setw(hn_weight_width) << hn_weight_stats.sd;
}

void printCommunityStats(const Statistic& community_size_stats,
                         const Statistic& community_pins_stats,
                         const Statistic& community_degree_stats) {
  // default double precision is 7
  const uint8_t double_width = 7;
  const uint8_t community_size_width = std::max(kahypar::math::digits(community_size_stats.max), double_width) + 4;
  const uint8_t community_pins_width = std::max(kahypar::math::digits(community_pins_stats.max), double_width) + 4;
  const uint8_t community_degree_width = std::max(kahypar::math::digits(community_degree_stats.max), double_width) + 4;

  LOG << "# HNs" << std::right << std::setw(community_size_width + 9)
      << "# Pins" << std::right << std::setw(community_pins_width + 8)
      << "Volume" << std::right << std::setw(community_degree_width + 8);
  LOG << "| min=" << std::left << std::setw(community_size_width) << community_size_stats.min
      << " | min=" << std::left << std::setw(community_pins_width) << community_pins_stats.min
      << " | min=" << std::left << std::setw(community_degree_width) << community_degree_stats.min;
  LOG << "| Q1 =" << std::left << std::setw(community_size_width) << community_size_stats.q1
      << " | Q1 =" << std::left << std::setw(community_pins_width) << community_pins_stats.q1
      << " | Q1 =" << std::left << std::setw(community_degree_width) << community_degree_stats.q1;
  LOG << "| med=" << std::left << std::setw(community_size_width) << community_size_stats.med
      << " | med=" << std::left << std::setw(community_pins_width) << community_pins_stats.med
      << " | med=" << std::left << std::setw(community_degree_width) << community_degree_stats.med;
  LOG << "| Q3 =" << std::left << std::setw(community_size_width) << community_size_stats.q3
      << " | Q3 =" << std::left << std::setw(community_pins_width) << community_pins_stats.q3
      << " | Q3 =" << std::left << std::setw(community_degree_width) << community_degree_stats.q3;
  LOG << "| max=" << std::left << std::setw(community_size_width) << community_size_stats.max
      << " | max=" << std::left << std::setw(community_pins_width) << community_pins_stats.max
      << " | max=" << std::left << std::setw(community_degree_width) << community_degree_stats.max;
  LOG << "| avg=" << std::left << std::setw(community_size_width) << community_size_stats.avg
      << " | avg=" << std::left << std::setw(community_pins_width) << community_pins_stats.avg
      << " | avg=" << std::left << std::setw(community_degree_width) << community_degree_stats.avg;
  LOG << "| sd =" << std::left << std::setw(community_size_width) << community_size_stats.sd
      << " | sd =" << std::left << std::setw(community_pins_width) << community_pins_stats.sd
      << " | sd =" << std::left << std::setw(community_degree_width) << community_degree_stats.sd;
}
}  // namespace internal

static inline void printBanner(const Context& context) {
  if (context.partition.verbose_output) {
    LOG << R"(+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++)";
    LOG << R"(+         __  __ _______       _  __     _    _       _____                   +)";
    LOG << R"(+        |  \/  |__   __|     | |/ /    | |  | |     |  __ \                  +)";
    LOG << R"(+        | \  / |  | |  ____  | ' / __ _| |__| |_   _| |__) |_ _ _ __         +)";
    LOG << R"(+        | |\/| |  | | |____| |  < / _` |  __  | | | |  ___/ _` | '__|        +)";
    LOG << R"(+        | |  | |  | |        | . \ (_| | |  | | |_| | |  | (_| | |           +)";
    LOG << R"(+        |_|  |_|  |_|        |_|\_\__,_|_|  |_|\__, |_|   \__,_|_|           +)";
    LOG << R"(+                                                __/ |                        +)";
    LOG << R"(+                                               |___/                         +)";
    LOG << R"(+          Karlsruhe Shared Memory Hypergraph Partitioning Framework          +)";
    LOG << R"(+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++)";
  }
}

template<typename HyperGraph>
inline void printHypergraphInfo(const HyperGraph& hypergraph, const std::string& name) {
  std::vector<HypernodeID> he_sizes;
  std::vector<HyperedgeWeight> he_weights;
  std::vector<HyperedgeID> hn_degrees;
  std::vector<HypernodeWeight> hn_weights;
  he_sizes.reserve(hypergraph.initialNumEdges());
  he_weights.reserve(hypergraph.initialNumEdges());
  hn_degrees.reserve(hypergraph.initialNumNodes());
  hn_weights.reserve(hypergraph.initialNumNodes());

  HypernodeID num_hypernodes = 0;
  const double avg_hn_degree = metrics::avgHypernodeDegree(hypergraph);
  double stdev_hn_degree = 0.0;
  for (const auto& hn : hypergraph.nodes()) {
    ++num_hypernodes;
    hn_degrees.push_back(hypergraph.nodeDegree(hn));
    hn_weights.push_back(hypergraph.nodeWeight(hn));
    stdev_hn_degree += (hypergraph.nodeDegree(hn) - avg_hn_degree) *
                       (hypergraph.nodeDegree(hn) - avg_hn_degree);
  }
  stdev_hn_degree = std::sqrt(stdev_hn_degree / (num_hypernodes - 1));

  HyperedgeID num_hyperedges = 0;
  HypernodeID num_pins = 0;
  const double avg_he_size = metrics::avgHyperedgeDegree(hypergraph);
  double stdev_he_size = 0.0;
  for (const auto& he : hypergraph.edges()) {
    ++num_hyperedges;
    num_pins += hypergraph.edgeSize(he);
    he_sizes.push_back(hypergraph.edgeSize(he));
    he_weights.push_back(hypergraph.edgeWeight(he));
    stdev_he_size += (hypergraph.edgeSize(he) - avg_he_size) *
                     (hypergraph.edgeSize(he) - avg_he_size);
  }
  stdev_he_size = std::sqrt(stdev_he_size / (num_hyperedges - 1));

  std::sort(he_sizes.begin(), he_sizes.end());
  std::sort(he_weights.begin(), he_weights.end());
  std::sort(hn_degrees.begin(), hn_degrees.end());
  std::sort(hn_weights.begin(), hn_weights.end());

  const double avg_hn_weight = std::accumulate(hn_weights.begin(), hn_weights.end(), 0.0) /
                               static_cast<double>(hn_weights.size());
  const double avg_he_weight = std::accumulate(he_weights.begin(), he_weights.end(), 0.0) /
                               static_cast<double>(he_weights.size());

  double stdev_hn_weight = 0.0;
  for (const HypernodeWeight& hn_weight : hn_weights) {
    stdev_hn_weight += (hn_weight - avg_hn_weight) * (hn_weight - avg_hn_weight);
  }
  stdev_hn_weight = std::sqrt(stdev_hn_weight / (num_hypernodes - 1));

  double stdev_he_weight = 0.0;
  for (const HyperedgeWeight& he_weight : he_weights) {
    stdev_he_weight += (he_weight - avg_he_weight) * (he_weight - avg_he_weight);
  }
  stdev_he_weight = std::sqrt(stdev_he_weight / (num_hyperedges - 1));

  LOG << "Hypergraph Information";
  LOG << "Name :" << name;
  LOG << "# HNs :" << num_hypernodes
      << "# HEs :" << num_hyperedges
      << "# pins:" << num_pins;

  internal::printHypergraphStats(
    internal::createStats(he_sizes, avg_he_size, stdev_he_size),
    internal::createStats(he_weights, avg_he_weight, stdev_he_weight),
    internal::createStats(hn_degrees, avg_hn_degree, stdev_hn_degree),
    internal::createStats(hn_weights, avg_hn_weight, stdev_hn_weight));

  // Print Memory Consumption
  utils::MemoryTreeNode hypergraph_memory_consumption("Hypergraph", utils::OutputType::MEGABYTE);
  hypergraph.memoryConsumption(&hypergraph_memory_consumption);
  hypergraph_memory_consumption.finalize();
  LOG << "\nHypergraph Memory Consumption";
  LOG << hypergraph_memory_consumption;
}

inline void printCommunityInformation(const Hypergraph& hypergraph) {
  PartitionID num_communities = hypergraph.numCommunities();
  std::vector<HypernodeID> community_sizes;
  std::vector<HypernodeID> community_pins;
  std::vector<HyperedgeID> community_degrees;
  community_sizes.reserve(num_communities);
  community_pins.reserve(num_communities);
  community_degrees.reserve(num_communities);

  for (PartitionID i = 0; i < num_communities; ++i) {
    community_sizes.push_back(hypergraph.numCommunityHypernodes(i));
    community_pins.push_back(hypergraph.numCommunityPins(i));
    community_degrees.push_back(hypergraph.communityDegree(i));
  }

  std::sort(community_sizes.begin(), community_sizes.end());
  std::sort(community_pins.begin(), community_pins.end());
  std::sort(community_degrees.begin(), community_degrees.end());

  const double avg_community_size = std::accumulate(community_sizes.begin(), community_sizes.end(), 0.0) /
                                    static_cast<double>(community_sizes.size());
  const double avg_community_pins = std::accumulate(community_pins.begin(), community_pins.end(), 0.0) /
                                    static_cast<double>(community_pins.size());
  const double avg_community_degree = std::accumulate(community_degrees.begin(), community_degrees.end(), 0.0) /
                                      static_cast<double>(community_degrees.size());

  double stdev_community_size = 0.0;
  double stdev_community_pins = 0.0;
  double stdev_community_degree = 0.0;
  for (PartitionID i = 0; i < num_communities; ++i) {
    stdev_community_size += (community_sizes[i] - avg_community_size) * (community_sizes[i] - avg_community_size);
    stdev_community_pins += (community_pins[i] - avg_community_pins) * (community_pins[i] - avg_community_pins);
    stdev_community_degree += (community_degrees[i] - avg_community_degree) * (community_degrees[i] - avg_community_degree);
  }
  stdev_community_size = std::sqrt(stdev_community_size / static_cast<double>(num_communities - 1));
  stdev_community_pins = std::sqrt(stdev_community_pins / static_cast<double>(num_communities - 1));
  stdev_community_degree = std::sqrt(stdev_community_degree / static_cast<double>(num_communities - 1));

  LOG << "Community Information:";
  LOG << "# Communities:" << num_communities;
  internal::printCommunityStats(
    internal::createStats(community_sizes, avg_community_size, stdev_community_size),
    internal::createStats(community_pins, avg_community_pins, stdev_community_pins),
    internal::createStats(community_degrees, avg_community_degree, stdev_community_degree));
}

inline void printPartSizesAndWeights(const PartitionedHypergraph<>& hypergraph, const Context& context) {
  HypernodeID max_part_size = 0;
  for (PartitionID i = 0; i != hypergraph.k(); ++i) {
    max_part_size = std::max(max_part_size, hypergraph.partSize(i));
  }
  const uint8_t part_digits = kahypar::math::digits(max_part_size);
  const uint8_t k_digits = kahypar::math::digits(hypergraph.k());
  for (PartitionID i = 0; i != hypergraph.k(); ++i) {
    bool is_imbalanced = hypergraph.partWeight(i) > context.partition.max_part_weights[i];
    if ( is_imbalanced ) std::cout << RED;
    std::cout << "|block " << std::left  << std::setw(k_digits) << i
              << std::setw(1) << "| = "  << std::right << std::setw(part_digits) << hypergraph.partSize(i)
              << std::setw(1) << "  w( "  << std::right << std::setw(k_digits) << i
              << std::setw(1) << " ) = "  << std::right << std::setw(part_digits) << hypergraph.partWeight(i)
              << std::setw(1) << "  max( " << std::right << std::setw(k_digits) << i
              << std::setw(1) << " ) = "  << std::right << std::setw(part_digits) << context.partition.max_part_weights[i]
              << std::endl;
    if ( is_imbalanced ) std::cout << END;
  }
}

static inline void printPartitioningResults(const PartitionedHypergraph<>& hypergraph,
                                            const Context& context,
                                            const std::string& description) {
  if (context.partition.verbose_output) {
    LOG << description;
    LOG << context.partition.objective << "      ="
        << metrics::objective(hypergraph, context.partition.objective);
    LOG << "imbalance =" << metrics::imbalance(hypergraph, context);
    LOG << "Part sizes and weights:";
    io::printPartSizesAndWeights(hypergraph, context);
    LOG << "";
  }
}

static inline void printContext(const Context& context) {
  if (context.partition.verbose_output) {
    LOG << context;
  }
}

static inline void printInputInformation(const Context& context, const Hypergraph& hypergraph) {
  if (context.partition.verbose_output) {
    LOG << "\n********************************************************************************";
    LOG << "*                                    Input                                     *";
    LOG << "********************************************************************************";
    io::printHypergraphInfo(hypergraph, context.partition.graph_filename.substr(
                              context.partition.graph_filename.find_last_of('/') + 1));
  }
}

static inline void printTopLevelPreprocessingBanner(const Context& context) {
  if (context.partition.verbose_output) {
    LOG << "\n********************************************************************************";
    LOG << "*                              Preprocessing...                                *";
    LOG << "********************************************************************************";
  }
}

static inline void printCoarseningBanner(const Context& context) {
  if (context.partition.verbose_output) {
    LOG << "********************************************************************************";
    LOG << "*                                Coarsening...                                 *";
    LOG << "********************************************************************************";
  }
}

static inline void printInitialPartitioningBanner(const Context& context) {
  if (context.partition.verbose_output) {
    LOG << "\n********************************************************************************";
    LOG << "*                           Initial Partitioning...                            *";
    LOG << "********************************************************************************";
  }
}

static inline void printLocalSearchBanner(const Context& context) {
  if (context.partition.verbose_output) {
    LOG << "\n********************************************************************************";
    LOG << "*                               Local Search...                                *";
    LOG << "********************************************************************************";
  }
}

inline void printObjectives(const PartitionedHypergraph<>& hypergraph,
                            const Context& context,
                            const std::chrono::duration<double>& elapsed_seconds) {
  LOG << "Objectives:";
  LOG << " Hyperedge Cut  (minimize) =" << metrics::hyperedgeCut(hypergraph);
  LOG << " SOED           (minimize) =" << metrics::soed(hypergraph);
  LOG << " (k-1)          (minimize) =" << metrics::km1(hypergraph);
  LOG << " Imbalance                 =" << metrics::imbalance(hypergraph, context);
  LOG << " Partitioning Time         =" << elapsed_seconds.count() << "s";
}

inline void printPartitioningResults(const PartitionedHypergraph<>& hypergraph,
                                     const Context& context,
                                     const std::chrono::duration<double>& elapsed_seconds) {
  unused(hypergraph);
  if (context.partition.verbose_output) {
    LOG << "\n********************************************************************************";
    LOG << "*                             Partitioning Result                              *";
    LOG << "********************************************************************************";
    printObjectives(hypergraph, context, elapsed_seconds);

    LOG << "\nPartition sizes and weights: ";
    printPartSizesAndWeights(hypergraph, context);

    LOG << "\nTimings:";
    LOG << utils::Timer::instance(context.partition.detailed_timings);
  }
}

static inline void printStripe() {
  LOG << "--------------------------------------------------------------------------------";
}
}  // namespace io
}  // namespace mt_kahypar

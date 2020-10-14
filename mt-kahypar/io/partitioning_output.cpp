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

#include "partitioning_output.h"

#include <vector>

#include "tbb/blocked_range.h"
#include "tbb/parallel_invoke.h"
#include "tbb/parallel_sort.h"
#include "tbb/parallel_reduce.h"
#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/memory_tree.h"
#include "mt-kahypar/utils/timer.h"


#include "kahypar/utils/math.h"



namespace mt_kahypar::io {
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

    template<typename T>
    double parallel_stdev(const std::vector<T>& data, const double avg, const size_t n) {
      return std::sqrt(tbb::parallel_reduce(
              tbb::blocked_range<size_t>(0UL, data.size()), 0.0,
              [&](tbb::blocked_range<size_t>& range, double init) -> double {
                double tmp_stdev = init;
                for ( size_t i = range.begin(); i < range.end(); ++i ) {
                  tmp_stdev += (data[i] - avg) * (data[i] - avg);
                }
                return tmp_stdev;
              }, std::plus<double>()) / ( n- 1 ));
    }

    template<typename T>
    double parallel_avg(const std::vector<T>& data, const size_t n) {
      return tbb::parallel_reduce(
              tbb::blocked_range<size_t>(0UL, data.size()), 0.0,
              [&](tbb::blocked_range<size_t>& range, double init) -> double {
                double tmp_avg = init;
                for ( size_t i = range.begin(); i < range.end(); ++i ) {
                  tmp_avg += static_cast<double>(data[i]);
                }
                return tmp_avg;
              }, std::plus<double>()) / static_cast<double>(n);
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

    static inline double avgHyperedgeDegree(const Hypergraph& hypergraph) {
      return static_cast<double>(hypergraph.initialNumPins()) / hypergraph.initialNumEdges();
    }

    static inline double avgHypernodeDegree(const Hypergraph& hypergraph) {
      return static_cast<double>(hypergraph.initialNumPins()) / hypergraph.initialNumNodes();
    }

  }  // namespace internal

  void printHypergraphInfo(const Hypergraph& hypergraph,
                                  const std::string& name,
                                  const bool show_memory_consumption) {
    std::vector<HypernodeID> he_sizes;
    std::vector<HyperedgeWeight> he_weights;
    std::vector<HyperedgeID> hn_degrees;
    std::vector<HypernodeWeight> hn_weights;

    tbb::parallel_invoke([&] {
      he_sizes.resize(hypergraph.initialNumEdges());
    }, [&] {
      he_weights.resize(hypergraph.initialNumEdges());
    }, [&] {
      hn_degrees.resize(hypergraph.initialNumNodes());
    }, [&] {
      hn_weights.resize(hypergraph.initialNumNodes());
    });

    HypernodeID num_hypernodes = hypergraph.initialNumNodes();
    const double avg_hn_degree = internal::avgHypernodeDegree(hypergraph);
    hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
      hn_degrees[hn] = hypergraph.nodeDegree(hn);
      hn_weights[hn] = hypergraph.nodeWeight(hn);
    });
    const double avg_hn_weight = internal::parallel_avg(hn_weights, num_hypernodes);
    const double stdev_hn_degree = internal::parallel_stdev(hn_degrees, avg_hn_degree, num_hypernodes);
    const double stdev_hn_weight = internal::parallel_stdev(hn_weights, avg_hn_weight, num_hypernodes);

    HyperedgeID num_hyperedges = hypergraph.initialNumEdges();
    HypernodeID num_pins = hypergraph.initialNumPins();
    const double avg_he_size = internal::avgHyperedgeDegree(hypergraph);
    hypergraph.doParallelForAllEdges([&](const HyperedgeID& he) {
      he_sizes[he] = hypergraph.edgeSize(he);
      he_weights[he] = hypergraph.edgeWeight(he);
    });
    const double avg_he_weight = internal::parallel_avg(he_weights, num_hyperedges);
    const double stdev_he_size = internal::parallel_stdev(he_sizes, avg_he_size, num_hyperedges);
    const double stdev_he_weight = internal::parallel_stdev(he_weights, avg_he_weight, num_hyperedges);

    tbb::enumerable_thread_specific<size_t> graph_edge_count(0);
    hypergraph.doParallelForAllEdges([&](const HyperedgeID& he) {
      if (hypergraph.edgeSize(he) == 2) {
        graph_edge_count.local() += 1;
      }
    });

    tbb::parallel_invoke([&] {
      tbb::parallel_sort(he_sizes.begin(), he_sizes.end());
    }, [&] {
      tbb::parallel_sort(he_weights.begin(), he_weights.end());
    }, [&] {
      tbb::parallel_sort(hn_degrees.begin(), hn_degrees.end());
    }, [&] {
      tbb::parallel_sort(hn_weights.begin(), hn_weights.end());
    });

    LOG << "Hypergraph Information";
    LOG << "Name :" << name;
    LOG << "# HNs :" << num_hypernodes
        << "# HEs :" << num_hyperedges
        << "# pins:" << num_pins
        << "# graph edges:" << graph_edge_count.combine(std::plus<>());

    internal::printHypergraphStats(
            internal::createStats(he_sizes, avg_he_size, stdev_he_size),
            internal::createStats(he_weights, avg_he_weight, stdev_he_weight),
            internal::createStats(hn_degrees, avg_hn_degree, stdev_hn_degree),
            internal::createStats(hn_weights, avg_hn_weight, stdev_hn_weight));

    if ( show_memory_consumption ) {
      // Print Memory Consumption
      utils::MemoryTreeNode hypergraph_memory_consumption("Hypergraph", utils::OutputType::MEGABYTE);
      hypergraph.memoryConsumption(&hypergraph_memory_consumption);
      hypergraph_memory_consumption.finalize();
      LOG << "\nHypergraph Memory Consumption";
      LOG << hypergraph_memory_consumption;
    }
  }

  void printPartWeightsAndSizes(const PartitionedHypergraph& hypergraph, const Context& context) {
    vec<HypernodeID> part_sizes(hypergraph.k(), 0);
    for (HypernodeID u : hypergraph.nodes()) {
      part_sizes[hypergraph.partID(u)]++;
    }
    HypernodeWeight max_part_weight = 0;
    HypernodeID max_part_size = 0;
    for (PartitionID i = 0; i < hypergraph.k(); ++i) {
      max_part_weight = std::max(max_part_weight, hypergraph.partWeight(i));
      max_part_size = std::max(max_part_size, part_sizes[i]);
    }
    const uint8_t part_digits = kahypar::math::digits(max_part_weight);
    const uint8_t k_digits = kahypar::math::digits(hypergraph.k());
    for (PartitionID i = 0; i != hypergraph.k(); ++i) {
      bool is_imbalanced =
              hypergraph.partWeight(i) > context.partition.max_part_weights[i] || hypergraph.partWeight(i) == 0;
      if ( is_imbalanced ) std::cout << RED;
      std::cout << "|block " << std::left  << std::setw(k_digits) << i
                << std::setw(1) << "| = "  << std::right << std::setw(part_digits) << part_sizes[i]
                << std::setw(1) << "  w( "  << std::right << std::setw(k_digits) << i
                << std::setw(1) << " ) = "  << std::right << std::setw(part_digits) << hypergraph.partWeight(i)
                << std::setw(1) << "  max( " << std::right << std::setw(k_digits) << i
                << std::setw(1) << " ) = "  << std::right << std::setw(part_digits) << context.partition.max_part_weights[i]
                << std::endl;
      if ( is_imbalanced ) std::cout << END;
    }
  }

  void printPartitioningResults(const PartitionedHypergraph& hypergraph,
                                              const Context& context,
                                              const std::string& description) {
    if (context.partition.verbose_output) {
      LOG << description;
      LOG << context.partition.objective << "      ="
          << metrics::objective(hypergraph, context.partition.objective);
      LOG << "imbalance =" << metrics::imbalance(hypergraph, context);
      LOG << "Part sizes and weights:";
      io::printPartWeightsAndSizes(hypergraph, context);
      LOG << "";
    }
  }

  void printContext(const Context& context) {
    if (context.partition.verbose_output) {
      LOG << context;
    }
  }

  void printMemoryPoolConsumption(const Context& context) {
    if ( context.partition.verbose_output && context.partition.show_memory_consumption ) {
      utils::MemoryTreeNode memory_pool_consumption("Memory Pool", utils::OutputType::MEGABYTE);
      parallel::MemoryPool::instance().memory_consumption(&memory_pool_consumption);
      memory_pool_consumption.finalize();
      LOG << "\n Memory Pool Consumption:";
      LOG << memory_pool_consumption << "\n";
      parallel::MemoryPool::instance().explain_optimizations();
    }
  }

  void printInputInformation(const Context& context, const Hypergraph& hypergraph) {
    if (context.partition.verbose_output) {
      LOG << "\n********************************************************************************";
      LOG << "*                                    Input                                     *";
      LOG << "********************************************************************************";
      io::printHypergraphInfo(hypergraph, context.partition.graph_filename.substr(
              context.partition.graph_filename.find_last_of('/') + 1),
                              context.partition.show_memory_consumption);
    }
  }

  void printTopLevelPreprocessingBanner(const Context& context) {
    if (context.partition.verbose_output) {
      LOG << "\n********************************************************************************";
      LOG << "*                              Preprocessing...                                *";
      LOG << "********************************************************************************";
    }
  }

  void printCoarseningBanner(const Context& context) {
    if (context.partition.verbose_output) {
      LOG << "********************************************************************************";
      LOG << "*                                Coarsening...                                 *";
      LOG << "********************************************************************************";
    }
  }

  void printInitialPartitioningBanner(const Context& context) {
    if (context.partition.verbose_output) {
      LOG << "\n********************************************************************************";
      LOG << "*                           Initial Partitioning...                            *";
      LOG << "********************************************************************************";
    }
  }

  void printLocalSearchBanner(const Context& context) {
    if (context.partition.verbose_output) {
      LOG << "\n********************************************************************************";
      LOG << "*                               Local Search...                                *";
      LOG << "********************************************************************************";
    }
  }

  void printVCycleBanner(const Context& context, const size_t vcycle_num) {
    if (context.partition.verbose_output) {
      LOG << "\n********************************************************************************";
      std::cout << "*                                  V-Cycle  " << vcycle_num;
      if ( vcycle_num < 10 ) {
        std::cout << "                                  *\n";
      } else {
        std::cout << "                                 *\n";
      }
      LOG << "********************************************************************************";
    }
  }

  void printObjectives(const PartitionedHypergraph& hypergraph,
                              const Context& context,
                              const std::chrono::duration<double>& elapsed_seconds) {
    LOG << "Objectives:";
    LOG << " Hyperedge Cut  (minimize) =" << metrics::hyperedgeCut(hypergraph);
    LOG << " SOED           (minimize) =" << metrics::soed(hypergraph);
    LOG << " (k-1)          (minimize) =" << metrics::km1(hypergraph);
    LOG << " Imbalance                 =" << metrics::imbalance(hypergraph, context);
    LOG << " Partitioning Time         =" << elapsed_seconds.count() << "s";
  }

  void printCutMatrix(const PartitionedHypergraph& hypergraph) {
    const PartitionID k = hypergraph.k();

    using MCell = parallel::IntegralAtomicWrapper<HyperedgeWeight>;
    using MCol = std::vector<MCell>;
    std::vector<MCol> cut_matrix(k, MCol(k, MCell(0)));

    hypergraph.doParallelForAllEdges([&](const HyperedgeID& he) {
      if ( hypergraph.connectivity(he) > 1 ) {
        const HyperedgeWeight edge_weight = hypergraph.edgeWeight(he);
        for ( const PartitionID& block_1 : hypergraph.connectivitySet(he) ) {
          for ( const PartitionID& block_2 : hypergraph.connectivitySet(he) ) {
            if ( block_1 < block_2 ) {
              cut_matrix[block_1][block_2] += edge_weight;
            }
          }
        }
      }
    });

    HyperedgeWeight max_cut = 0;
    for ( PartitionID block_1 = 0; block_1 < k; ++block_1 ) {
      for ( PartitionID block_2 = block_1 + 1; block_2 < k; ++block_2 ) {
        max_cut = std::max(max_cut, cut_matrix[block_1][block_2].load());
      }
    }

    // HEADER
    const uint8_t column_width = std::max(kahypar::math::digits(max_cut) + 2, 5);
    std::cout << std::right << std::setw(column_width) << "Block";
    for ( PartitionID block = 0; block < k; ++block ) {
      std::cout << std::right << std::setw(column_width) << block;
    }
    std::cout << std::endl;

    // CUT MATRIX
    for ( PartitionID block_1 = 0; block_1 < k; ++block_1 ) {
      std::cout << std::right << std::setw(column_width) << block_1;
      for ( PartitionID block_2 = 0; block_2 < k; ++block_2 ) {
        std::cout << std::right << std::setw(column_width) << cut_matrix[block_1][block_2].load();
      }
      std::cout << std::endl;
    }
  }

  void printPotentialPositiveGainMoveMatrix(const PartitionedHypergraph& hypergraph) {
    const PartitionID k = hypergraph.k();

    using MCell = parallel::IntegralAtomicWrapper<HyperedgeWeight>;
    using MCol = std::vector<MCell>;
    std::vector<MCol> positive_gains(k, MCol(k, MCell(0)));

    tbb::enumerable_thread_specific<std::vector<Gain>> local_gain(k, 0);
    hypergraph.doParallelForAllNodes([&](const HypernodeID hn) {
      // Calculate gain to all blocks of the partition
      std::vector<Gain>& tmp_scores = local_gain.local();
      const PartitionID from = hypergraph.partID(hn);
      Gain internal_weight = 0;
      for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
        HypernodeID pin_count_in_from_part = hypergraph.pinCountInPart(he, from);
        HyperedgeWeight he_weight = hypergraph.edgeWeight(he);

        if ( pin_count_in_from_part > 1 ) {
          internal_weight += he_weight;
        }

        for (const PartitionID& to : hypergraph.connectivitySet(he)) {
          if (from != to) {
            tmp_scores[to] -= he_weight;
          }
        }
      }

      for (PartitionID to = 0; to < k; ++to) {
        if (from != to) {
          Gain score = tmp_scores[to] + internal_weight;
          if ( score < 0 ) {
            positive_gains[from][to] += std::abs(score);
          }
        }
        tmp_scores[to] = 0;
      }
    });


    HyperedgeWeight max_gain = 0;
    for ( PartitionID block_1 = 0; block_1 < k; ++block_1 ) {
      for ( PartitionID block_2 = block_1 + 1; block_2 < k; ++block_2 ) {
        max_gain = std::max(max_gain, positive_gains[block_1][block_2].load());
      }
    }

    // HEADER
    const uint8_t column_width = std::max(kahypar::math::digits(max_gain) + 2, 5);
    std::cout << std::right << std::setw(column_width) << "Block";
    for ( PartitionID block = 0; block < k; ++block ) {
      std::cout << std::right << std::setw(column_width) << block;
    }
    std::cout << std::endl;

    // CUT MATRIX
    for ( PartitionID block_1 = 0; block_1 < k; ++block_1 ) {
      std::cout << std::right << std::setw(column_width) << block_1;
      for ( PartitionID block_2 = 0; block_2 < k; ++block_2 ) {
        std::cout << std::right << std::setw(column_width) << positive_gains[block_1][block_2].load();
      }
      std::cout << std::endl;
    }
  }

  void printConnectedCutHyperedgeAnalysis(const PartitionedHypergraph& hypergraph) {
    std::vector<bool> visited_he(hypergraph.initialNumEdges(), false);
    std::vector<HyperedgeWeight> connected_cut_hyperedges;

    auto analyse_component = [&](const HyperedgeID he) {
      HyperedgeWeight component_weight = 0;
      std::vector<HyperedgeID> s;
      s.push_back(he);
      visited_he[he] = true;

      while ( !s.empty() ) {
        const HyperedgeID e = s.back();
        s.pop_back();
        component_weight += hypergraph.edgeWeight(e);

        for ( const HypernodeID& pin : hypergraph.pins(e) ) {
          for ( const HyperedgeID& tmp_e : hypergraph.incidentEdges(pin) ) {
            if ( !visited_he[tmp_e] && hypergraph.connectivity(tmp_e) > 1 ) {
              s.push_back(tmp_e);
              visited_he[tmp_e] = true;
            }
          }
        }
      }

      return component_weight;
    };

    for ( const HyperedgeID& he : hypergraph.edges() ) {
      if ( hypergraph.connectivity(he) > 1 && !visited_he[he] ) {
        connected_cut_hyperedges.push_back(analyse_component(he));
      }
    }
    std::sort(connected_cut_hyperedges.begin(), connected_cut_hyperedges.end());

    LOG << "Num Connected Cut Hyperedges =" << connected_cut_hyperedges.size();
    LOG << "Min Component Weight         =" << connected_cut_hyperedges[0];
    LOG << "Median Component Weight      =" << connected_cut_hyperedges[connected_cut_hyperedges.size() / 2];
    LOG << "Max Component Weight         =" << connected_cut_hyperedges.back();
    LOG << "Component Weight Vector:";
    std::cout << "(";
    for ( const HyperedgeWeight& weight : connected_cut_hyperedges ) {
      std::cout << weight << ",";
    }
    std::cout << "\b)" << std::endl;
  }

  void printPartitioningResults(const PartitionedHypergraph& hypergraph,
                                       const Context& context,
                                       const std::chrono::duration<double>& elapsed_seconds) {
    if (context.partition.verbose_output) {
      LOG << "\n********************************************************************************";
      LOG << "*                             Partitioning Result                              *";
      LOG << "********************************************************************************";
      printObjectives(hypergraph, context, elapsed_seconds);

      if ( context.partition.show_advanced_cut_analysis ) {
        LOG << "\nCut Matrix: ";
        printCutMatrix(hypergraph);

        LOG << "\nPotential Positive Gain Move Matrix: ";
        printPotentialPositiveGainMoveMatrix(hypergraph);

        LOG << "\nConnected Cut Hyperedge Analysis: ";
        printConnectedCutHyperedgeAnalysis(hypergraph);
      }

      LOG << "\nPartition sizes and weights: ";
      printPartWeightsAndSizes(hypergraph, context);

      LOG << "\nTimings:";
      utils::Timer::instance().setMaximumOutputDepth(context.partition.timings_output_depth);
      LOG << utils::Timer::instance(context.partition.show_detailed_timings);
    }
  }

  void printStripe() {
    LOG << "--------------------------------------------------------------------------------";
  }

  void printBanner() {
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

  namespace internal {
    void printCommunityStats(const Statistic& community_size_stats,
                             const Statistic& community_pins_stats,
                             const Statistic& community_degree_stats) {
      // default double precision is 7
      const uint8_t double_width = 7;
      const uint8_t community_size_width = std::max(kahypar::math::digits(community_size_stats.max), double_width) + 4;
      const uint8_t community_pins_width = std::max(kahypar::math::digits(community_pins_stats.max), double_width) + 4;
      const uint8_t community_degree_width = std::max(kahypar::math::digits(community_degree_stats.max), double_width) + 4;

      LOG << "# HNs Per Community" << std::right << std::setw(community_size_width + 2)
          << "# Internal Pins" << std::right << std::setw(community_pins_width + 8 + 4)
          << "Internal Degree Sum" << std::right << std::setw(community_degree_width + 8);
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
  }



  void printCommunityInformation(const Hypergraph& hypergraph) {

    PartitionID num_communities =
            tbb::parallel_reduce(
                    tbb::blocked_range<HypernodeID>(ID(0), hypergraph.initialNumNodes()),
                    0, [&](const tbb::blocked_range<HypernodeID>& range, PartitionID init) {
              PartitionID my_range_num_communities = init;
              for (HypernodeID hn = range.begin(); hn < range.end(); ++hn) {
                if ( hypergraph.nodeIsEnabled(hn) ) {
                  my_range_num_communities = std::max(my_range_num_communities, hypergraph.communityID(hn) + 1);
                }
              }
              return my_range_num_communities;
            },
            [](const PartitionID lhs, const PartitionID rhs) {
              return std::max(lhs, rhs);
            });
    num_communities = std::max(num_communities, 1);

    std::vector<size_t> nodes_per_community(num_communities, 0);
    std::vector<size_t> internal_pins(num_communities, 0);
    std::vector<size_t> internal_degree(num_communities, 0);

    auto reduce_nodes = [&] {
      tbb::enumerable_thread_specific< vec< std::pair<size_t, size_t> > > ets_nodes(num_communities, std::make_pair(0UL, 0UL));
      hypergraph.doParallelForAllNodes([&](const HypernodeID u) {
        const PartitionID cu = hypergraph.communityID(u);
        ets_nodes.local()[cu].first++;
        ets_nodes.local()[cu].second += hypergraph.nodeDegree(u);
      });

      for (const auto& x : ets_nodes) {
        for (PartitionID i = 0; i < num_communities; ++i) {
          nodes_per_community[i] += x[i].first;
          internal_degree[i] += x[i].second;
        }
      }
    };

    auto reduce_hyperedges = [&] {
      tbb::enumerable_thread_specific< vec<size_t> > ets_pins(num_communities, 0);
      hypergraph.doParallelForAllEdges([&](const HyperedgeID he) {
        auto& pin_counter = ets_pins.local();
        for (const HypernodeID pin : hypergraph.pins(he)) {
          pin_counter[ hypergraph.communityID(pin) ]++;
        }
      });

      for (const auto& x : ets_pins) {
        for (PartitionID i = 0; i < num_communities; ++i) {
          internal_pins[i] += x[i];
        }
      }
    };

    tbb::parallel_invoke(reduce_nodes, reduce_hyperedges);

    std::sort(nodes_per_community.begin(), nodes_per_community.end());
    std::sort(internal_pins.begin(), internal_pins.end());
    std::sort(internal_degree.begin(), internal_degree.end());

    auto square = [&](size_t x) { return x * x; };

    auto avg_and_std_dev = [&](const std::vector<size_t>& v) {
      const double avg = std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(v.size());
      double std_dev = 0.0;
      for (size_t x : v) {
        std_dev += square(x - avg);
      }
      std_dev = std::sqrt(std_dev / static_cast<double>(v.size() - 1));
      return std::make_pair(avg, std_dev);
    };

    auto [avg_nodes, std_dev_nodes] = avg_and_std_dev(nodes_per_community);
    auto [avg_pins, std_dev_pins] = avg_and_std_dev(internal_pins);
    auto [avg_deg, std_dev_deg] = avg_and_std_dev(internal_degree);

    LOG << "# Communities :" << num_communities;

    internal::printCommunityStats(
            internal::createStats(nodes_per_community, avg_nodes, std_dev_nodes),
            internal::createStats(internal_pins, avg_pins, std_dev_pins),
            internal::createStats(internal_degree, avg_deg, std_dev_deg)
            );
  }


} // namespace mt_kahypar::io
/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@student.kit.edu>
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#pragma once

#include "tbb/parallel_sort.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/utils/hypergraph_statistics.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

template<typename TypeTraits>
class Level {

  using Hypergraph = typename TypeTraits::Hypergraph;

public:
  explicit Level(Hypergraph&& contracted_hypergraph,
                 parallel::scalable_vector<HypernodeID>&& communities,
                 double coarsening_time) :
    _contracted_hypergraph(std::move(contracted_hypergraph)),
    _communities(std::move(communities)),
    _coarsening_time(coarsening_time) { }

  Hypergraph& contractedHypergraph() {
    return _contracted_hypergraph;
  }

  const Hypergraph& contractedHypergraph() const {
    return _contracted_hypergraph;
  }

  // ! Maps a global vertex id of the representative hypergraph
  // ! to its global vertex id in the contracted hypergraph
  HypernodeID mapToContractedHypergraph(const HypernodeID hn) const {
    ASSERT(hn < _communities.size());
    return _communities[hn];
  }

  double coarseningTime() const {
    return _coarsening_time;
  }

  void freeInternalData() {
    tbb::parallel_invoke([&] {
      _contracted_hypergraph.freeInternalData();
    }, [&] {
      parallel::free(_communities);
    });
  }

private:
  // ! Contracted Hypergraph
  Hypergraph _contracted_hypergraph;
  // ! Defines the communities that are contracted
  // ! in the coarse hypergraph
  parallel::scalable_vector<HypernodeID> _communities;
  // ! Time to create the coarsened hypergraph
  // ! (includes coarsening + contraction time)
  double _coarsening_time;
};

template<typename TypeTraits>
class UncoarseningData {

  using Hypergraph = typename TypeTraits::Hypergraph;
  using HypergraphFactory = typename Hypergraph::Factory;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
  using ParallelHyperedge = typename Hypergraph::ParallelHyperedge;

public:
  explicit UncoarseningData(bool n_level, Hypergraph& hg, const Context& context) :
    nlevel(n_level),
    _hg(hg),
    _context(context) {
      if (n_level) {
        compactified_hg = std::make_unique<Hypergraph>();
        compactified_phg = std::make_unique<PartitionedHypergraph>();
      } else {
        size_t estimated_number_of_levels = UL(1);
        if ( hg.initialNumNodes() > context.coarsening.contraction_limit ) {
          estimated_number_of_levels = std::ceil( std::log2(
              static_cast<double>(hg.initialNumNodes()) /
              static_cast<double>(context.coarsening.contraction_limit)) /
            std::log2(context.coarsening.maximum_shrink_factor) ) + UL(1);
        }
        hierarchy.reserve(estimated_number_of_levels);
      }
      is_phg_initialized = false;
      partitioned_hg = std::make_unique<PartitionedHypergraph>();
    }

  ~UncoarseningData() noexcept {
    tbb::parallel_for(UL(0), hierarchy.size(), [&](const size_t i) {
      (hierarchy)[i].freeInternalData();
    }, tbb::static_partitioner());
  }

  void setPartitionedHypergraph(PartitionedHypergraph&& phg) {
    ASSERT(!is_phg_initialized);
    partitioned_hg = std::make_unique<PartitionedHypergraph>(std::move(phg));
    is_phg_initialized = true;
  }

  void finalizeCoarsening() {
    utils::Timer& timer = utils::Utilities::instance().getTimer(_context.utility_id);
    if (nlevel) {
      // Create compactified hypergraph containing only enabled vertices and hyperedges
      // with consecutive IDs => Less complexity in initial partitioning.
      timer.start_timer("compactify_hypergraph", "Compactify Hypergraph");
      auto compactification = HypergraphFactory::compactify(_hg);
      *compactified_hg = std::move(compactification.first);
      compactified_hn_mapping = std::move(compactification.second);
      *compactified_phg = PartitionedHypergraph(_context.partition.k, *compactified_hg, parallel_tag_t());
      timer.stop_timer("compactify_hypergraph");
    } else {
      timer.start_timer("finalize_multilevel_hierarchy", "Finalize Multilevel Hierarchy");
      // Free memory of temporary contraction buffer and
      // release coarsening memory in memory pool
      if (!hierarchy.empty()) {
        hierarchy.back().contractedHypergraph().freeTmpContractionBuffer();
      } else {
        _hg.freeTmpContractionBuffer();
      }
      if (_context.type == ContextType::main) {
        parallel::MemoryPool::instance().release_mem_group("Coarsening");
      }

      // Construct partitioned hypergraph for initial partitioning
      if ( !is_phg_initialized ) {
        *partitioned_hg = PartitionedHypergraph(_context.partition.k, _hg, parallel_tag_t());
      }
      if (!hierarchy.empty()) {
        partitioned_hg->setHypergraph(hierarchy.back().contractedHypergraph());
      }
      is_phg_initialized = true;
      timer.stop_timer("finalize_multilevel_hierarchy");
    }
    is_finalized = true;
  }

  void performMultilevelContraction(
          parallel::scalable_vector<HypernodeID>&& communities, bool deterministic,
          const HighResClockTimepoint& round_start) {

    auto output_stats = [](const Hypergraph& hypergraph) {
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
      const double avg_hn_degree = utils::avgHypernodeDegree(hypergraph);
      hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
        hn_degrees[hn] = hypergraph.nodeDegree(hn);
        hn_weights[hn] = hypergraph.nodeWeight(hn);
      });
      const double avg_hn_weight = utils::parallel_avg(hn_weights, num_hypernodes);
      const double stdev_hn_degree = utils::parallel_stdev(hn_degrees, avg_hn_degree, num_hypernodes);
      const double stdev_hn_weight = utils::parallel_stdev(hn_weights, avg_hn_weight, num_hypernodes);

      HyperedgeID num_hyperedges = hypergraph.initialNumEdges();
      HypernodeID num_pins = hypergraph.initialNumPins();
      const double avg_he_size = utils::avgHyperedgeDegree(hypergraph);
      hypergraph.doParallelForAllEdges([&](const HyperedgeID& he) {
        he_sizes[he] = hypergraph.edgeSize(he);
        he_weights[he] = hypergraph.edgeWeight(he);
      });
      const double avg_he_weight = utils::parallel_avg(he_weights, num_hyperedges);
      const double stdev_he_size = utils::parallel_stdev(he_sizes, avg_he_size, num_hyperedges);
      const double stdev_he_weight = utils::parallel_stdev(he_weights, avg_he_weight, num_hyperedges);

      tbb::parallel_invoke([&] {
        tbb::parallel_sort(he_sizes.begin(), he_sizes.end());
      }, [&] {
        tbb::parallel_sort(he_weights.begin(), he_weights.end());
      }, [&] {
        tbb::parallel_sort(hn_degrees.begin(), hn_degrees.end());
      }, [&] {
        tbb::parallel_sort(hn_weights.begin(), hn_weights.end());
      });

      auto print_histogramm = [](const auto& sorted_values) {
        uint32_t upper = 1;
        size_t i = 0;
        while (i < sorted_values.size()) {
          size_t count = 0;
          while (sorted_values[i] <= upper) {
            ++count;
            ++i;
          }
          std::cout << count << " [<=" << upper << "] ";
          upper *= 2;
        }
        LOG << "";
      };

      LOG << "";
      LOG << V(num_hypernodes) << V(num_hyperedges);
      LOG << "HE weights distribution:";
      print_histogramm(he_weights);
      LOG << "HN degrees distribution:";
      print_histogramm(hn_degrees);
      LOG << "HN weights distribution:";
      print_histogramm(hn_weights);
    };

    ASSERT(!is_finalized);
    Hypergraph& current_hg = hierarchy.empty() ? _hg : hierarchy.back().contractedHypergraph();
    ASSERT(current_hg.initialNumNodes() == communities.size());
    if (hierarchy.empty()) {
      output_stats(current_hg);
    }
    Hypergraph contracted_hg = current_hg.contract(communities, deterministic);
    output_stats(contracted_hg);

    const HighResClockTimepoint round_end = std::chrono::high_resolution_clock::now();
    const double elapsed_time = std::chrono::duration<double>(round_end - round_start).count();
    hierarchy.emplace_back(std::move(contracted_hg), std::move(communities), elapsed_time);
  }

  PartitionedHypergraph& coarsestPartitionedHypergraph() {
    if (nlevel) {
      return *compactified_phg;
    } else {
      return *partitioned_hg;
    }
  }

  // Multilevel Data
  vec<Level<TypeTraits>> hierarchy;

  // NLevel Data
  // ! Once coarsening terminates we generate a compactified hypergraph
  // ! containing only enabled vertices and hyperedges within a consecutive
  // ! ID range, which is then used for initial partitioning
  std::unique_ptr<Hypergraph> compactified_hg;
  // ! Mapping from vertex IDs of the original hypergraph to the IDs
  // ! in the compactified hypergraph
  vec<HypernodeID> compactified_hn_mapping;
  // ! Compactified partitioned hypergraph
  std::unique_ptr<PartitionedHypergraph> compactified_phg;
  // ! Contains timings how long a coarsening pass takes for each round
  vec<vec<ParallelHyperedge>> removed_hyperedges_batches;
  // ! Removed single-pin and parallel nets.
  // ! All hyperedges that are contained in one vector must be restored once
  // ! we completly processed a vector of batches.
  vec<double> round_coarsening_times;

  // Both
  bool is_phg_initialized;
  std::unique_ptr<PartitionedHypergraph> partitioned_hg;
  bool is_finalized = false;
  bool nlevel;

private:
  Hypergraph& _hg;
  const Context& _context;
};

typedef struct uncoarsening_data_s uncoarsening_data_t;

namespace uncoarsening {
  template<typename TypeTraits>
  uncoarsening_data_t* to_pointer(UncoarseningData<TypeTraits>& ip_data) {
    return reinterpret_cast<uncoarsening_data_t*>(&ip_data);
  }

  template<typename TypeTraits>
  UncoarseningData<TypeTraits>& to_reference(uncoarsening_data_t* ptr) {
    return *reinterpret_cast<UncoarseningData<TypeTraits>*>(ptr);
  }
}

}

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

#include "tbb/task_group.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/rebalancing/rebalancer.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/utils/progress_bar.h"
#include "mt-kahypar/utils/stats.h"

namespace mt_kahypar {

class MultilevelCoarsenerBase {
 private:

  static constexpr bool debug = false;

  class Level {

   public:
    explicit Level(Hypergraph&& contracted_hypergraph,
                   parallel::scalable_vector<HypernodeID>&& communities,
                   double coarsening_time) :
      _representative_hypergraph(nullptr),
      _contracted_hypergraph(std::move(contracted_hypergraph)),
      _contracted_partitioned_hypergraph(),
      _communities(std::move(communities)),
      _coarsening_time(coarsening_time) { }

    void setRepresentativeHypergraph(PartitionedHypergraph* representative_hypergraph) {
      _representative_hypergraph = representative_hypergraph;
    }

    PartitionedHypergraph& representativeHypergraph() {
      ASSERT(_representative_hypergraph);
      return *_representative_hypergraph;
    }

    Hypergraph& contractedHypergraph() {
      return _contracted_hypergraph;
    }

    PartitionedHypergraph& contractedPartitionedHypergraph() {
      return _contracted_partitioned_hypergraph;
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
        _contracted_partitioned_hypergraph.freeInternalData();
      }, [&] {
        parallel::free(_communities);
      });
    }

   private:
    // ! Hypergraph on the next finer level
    PartitionedHypergraph* _representative_hypergraph;
    // ! Contracted Hypergraph
    Hypergraph _contracted_hypergraph;
    // ! Partitioned Hypergraph
    PartitionedHypergraph _contracted_partitioned_hypergraph;
    // ! Defines the communities that are contracted
    // ! in the coarse hypergraph
    parallel::scalable_vector<HypernodeID> _communities;
    // ! Time to create the coarsened hypergraph
    // ! (includes coarsening + contraction time)
    double _coarsening_time;
  };

 public:
  MultilevelCoarsenerBase(Hypergraph& hypergraph,
                          const Context& context,
                          const TaskGroupID task_group_id,
                          const bool top_level) :
          _is_finalized(false),
          _hg(hypergraph),
          _partitioned_hg(),
          _context(context),
          _task_group_id(task_group_id),
          _top_level(top_level),
          _hierarchy() {
    size_t estimated_number_of_levels = 1UL;
    if ( _hg.initialNumNodes() > _context.coarsening.contraction_limit ) {
      estimated_number_of_levels = std::ceil( std::log2(
        static_cast<double>(_hg.initialNumNodes()) /
        static_cast<double>(_context.coarsening.contraction_limit)) /
        std::log2(_context.coarsening.maximum_shrink_factor) ) + 1UL;
    }
    _hierarchy.reserve(estimated_number_of_levels);
  }

  MultilevelCoarsenerBase(const MultilevelCoarsenerBase&) = delete;
  MultilevelCoarsenerBase(MultilevelCoarsenerBase&&) = delete;
  MultilevelCoarsenerBase & operator= (const MultilevelCoarsenerBase &) = delete;
  MultilevelCoarsenerBase & operator= (MultilevelCoarsenerBase &&) = delete;

  virtual ~MultilevelCoarsenerBase() throw () {
    tbb::parallel_for(0UL, _hierarchy.size(), [&](const size_t i) {
      _hierarchy[i].freeInternalData();
    }, tbb::static_partitioner());
  }

 protected:

  HypernodeID currentNumNodes() const {
    if ( _hierarchy.empty() ) {
      return _hg.initialNumNodes();
    } else {
      return _hierarchy.back().contractedHypergraph().initialNumNodes();
    }
  }

  Hypergraph& currentHypergraph() {
    if ( _hierarchy.empty() ) {
      return _hg;
    } else {
      return _hierarchy.back().contractedHypergraph();
    }
  }

  PartitionedHypergraph& currentPartitionedHypergraph() {
    ASSERT(_is_finalized);
    if ( _hierarchy.empty() ) {
      return _partitioned_hg;
    } else {
      return _hierarchy.back().contractedPartitionedHypergraph();
    }
  }

  void finalize();

  void performMultilevelContraction(parallel::scalable_vector<HypernodeID>&& communities, const HighResClockTimepoint& round_start);

  PartitionedHypergraph&& doUncoarsen(std::unique_ptr<IRefiner>& label_propagation, std::unique_ptr<IRefiner>& fm);

 protected:

  kahypar::Metrics computeMetrics(PartitionedHypergraph& phg) {
    HyperedgeWeight cut = 0;
    HyperedgeWeight km1 = 0;
    tbb::parallel_invoke([&] {
      cut = metrics::hyperedgeCut(phg);
    }, [&] {
      km1 = metrics::km1(phg);
    });
    return { cut, km1,  metrics::imbalance(phg, _context) };
  }

  kahypar::Metrics initialize(PartitionedHypergraph& current_hg) {
    kahypar::Metrics current_metrics = computeMetrics(current_hg);
    int64_t num_nodes = current_hg.initialNumNodes();
    int64_t num_edges = current_hg.initialNumEdges();
    utils::Stats::instance().add_stat("initial_num_nodes", num_nodes);
    utils::Stats::instance().add_stat("initial_num_edges", num_edges);
    utils::Stats::instance().add_stat("initial_cut", current_metrics.cut);
    utils::Stats::instance().add_stat("initial_km1", current_metrics.km1);
    utils::Stats::instance().add_stat("initial_imbalance", current_metrics.imbalance);
    return current_metrics;
  }

  void refine(PartitionedHypergraph& partitioned_hypergraph,
              std::unique_ptr<IRefiner>& label_propagation,
              std::unique_ptr<IRefiner>& fm,
              kahypar::Metrics& current_metrics,
              const double time_limit);

  double refinementTimeLimit(const Level& hierarchy) const {
    if ( _context.refinement.fm.time_limit_factor != std::numeric_limits<double>::max() ) {
      const double time_limit_factor = std::max(1.0,  _context.refinement.fm.time_limit_factor * _context.partition.k);
      return std::max(5.0, time_limit_factor * hierarchy.coarseningTime());
    } else {
      return std::numeric_limits<double>::max();
    }
  }

  bool _is_finalized;
  Hypergraph& _hg;
  PartitionedHypergraph _partitioned_hg;
  const Context& _context;
  const TaskGroupID _task_group_id;
  const bool _top_level;
  vec<Level> _hierarchy;
};
}  // namespace mt_kahypar

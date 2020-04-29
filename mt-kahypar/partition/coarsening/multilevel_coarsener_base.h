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
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/rebalancing/rebalancer.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/utils/progress_bar.h"
#include "mt-kahypar/utils/stats.h"

#include "mt-kahypar/partition/refinement/fm/multitry_kway_fm.h"  // for now. include it in registries later

namespace mt_kahypar {

class MultilevelCoarsenerBase {
 private:

  static constexpr bool debug = true;

  class Hierarchy {

   public:
    explicit Hierarchy(Hypergraph&& contracted_hypergraph,
                       parallel::scalable_vector<HypernodeID>&& communities) :
      _representative_hypergraph(nullptr),
      _contracted_hypergraph(std::move(contracted_hypergraph)),
      _contracted_partitioned_hypergraph(),
      _communities(std::move(communities)) { }

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
    _hierarchies() {
    size_t estimated_number_of_levels = 1UL;
    if ( _hg.initialNumNodes() > _context.coarsening.contraction_limit ) {
      estimated_number_of_levels = std::ceil( std::log2(
        static_cast<double>(_hg.initialNumNodes()) /
        static_cast<double>(_context.coarsening.contraction_limit)) /
        std::log2(_context.coarsening.maximum_shrink_factor) ) + 1UL;
    }
    _hierarchies.reserve(estimated_number_of_levels);
  }

  MultilevelCoarsenerBase(const MultilevelCoarsenerBase&) = delete;
  MultilevelCoarsenerBase(MultilevelCoarsenerBase&&) = delete;
  MultilevelCoarsenerBase & operator= (const MultilevelCoarsenerBase &) = delete;
  MultilevelCoarsenerBase & operator= (MultilevelCoarsenerBase &&) = delete;

  virtual ~MultilevelCoarsenerBase() throw () {
    tbb::parallel_for(0UL, _hierarchies.size(), [&](const size_t i) {
      _hierarchies[i].freeInternalData();
    }, tbb::static_partitioner());
  }

 protected:

  HypernodeID currentNumNodes() const {
    if ( _hierarchies.empty() ) {
      return _hg.initialNumNodes();
    } else {
      return _hierarchies.back().contractedHypergraph().initialNumNodes();
    }
  }

  Hypergraph& currentHypergraph() {
    if ( _hierarchies.empty() ) {
      return _hg;
    } else {
      return _hierarchies.back().contractedHypergraph();
    }
  }

  PartitionedHypergraph& currentPartitionedHypergraph() {
    ASSERT(_is_finalized);
    if ( _hierarchies.empty() ) {
      return _partitioned_hg;
    } else {
      return _hierarchies.back().contractedPartitionedHypergraph();
    }
  }

  void finalize() {
    utils::Timer::instance().start_timer("finalize_multilevel_hierarchy", "Finalize Multilevel Hierarchy");
    // Free memory of temporary contraction buffer and
    // release coarsening memory in memory pool
    currentHypergraph().freeTmpContractionBuffer();
    if ( _top_level ) {
      parallel::MemoryPool::instance().release_mem_group("Coarsening");
    }

    // Construct top level partitioned hypergraph (memory is taken from memory pool)
    _partitioned_hg = PartitionedHypergraph(
      _context.partition.k, _task_group_id, _hg);

    // Construct partitioned hypergraphs parallel
    tbb::task_group group;
    // Construct partitioned hypergraph for each coarsened hypergraph in the hierarchy
    for ( size_t i = 0; i < _hierarchies.size(); ++i ) {
      group.run([&, i] {
        _hierarchies[i].contractedPartitionedHypergraph() = PartitionedHypergraph(
           _context.partition.k, _task_group_id, _hierarchies[i].contractedHypergraph());
      });
    }
    group.wait();

    // Set the representative partitioned hypergraph for each hypergraph
    // in the hierarchy
    if ( _hierarchies.size() > 0 ) {
      _hierarchies[0].setRepresentativeHypergraph(&_partitioned_hg);
      for ( size_t i = 1; i < _hierarchies.size(); ++i ) {
        _hierarchies[i].setRepresentativeHypergraph(&_hierarchies[i - 1].contractedPartitionedHypergraph());
      }
    }
    _is_finalized = true;
    utils::Timer::instance().stop_timer("finalize_multilevel_hierarchy");
  }

  void performMultilevelContraction(parallel::scalable_vector<HypernodeID>&& communities) {
    ASSERT(!_is_finalized);
    Hypergraph& current_hg = currentHypergraph();
    ASSERT(current_hg.initialNumNodes() == communities.size());
    Hypergraph contracted_hg = current_hg.contract(communities, _task_group_id);
    _hierarchies.emplace_back(std::move(contracted_hg), std::move(communities));
  }

  PartitionedHypergraph&& doUncoarsen(std::unique_ptr<IRefiner>& label_propagation) {
    PartitionedHypergraph& coarsest_hg = currentPartitionedHypergraph();
    kahypar::Metrics current_metrics = initialize(coarsest_hg);

    refinement::MultiTryKWayFM fm_refiner(_context, _task_group_id, _hg.initialNumNodes(), _hg.initialNumEdges());

    utils::ProgressBar uncontraction_progress(_hg.initialNumNodes(),
      _context.partition.objective == kahypar::Objective::km1 ? current_metrics.km1 : current_metrics.cut,
      _context.partition.verbose_output && _context.partition.enable_progress_bar);
    uncontraction_progress += coarsest_hg.initialNumNodes();

    // Refine Coarsest Partitioned Hypergraph
    LOG << "Refine coarsest HG" << V(_context.partition.k) << V(coarsest_hg.initialNumNodes());
    refine(coarsest_hg, label_propagation, fm_refiner, current_metrics);

    for ( int i = _hierarchies.size() - 1; i >= 0; --i ) {
      // Project partition to next level finer hypergraph
      utils::Timer::instance().start_timer("projecting_partition", "Projecting Partition");
      PartitionedHypergraph& representative_hg = _hierarchies[i].representativeHypergraph();
      PartitionedHypergraph& contracted_hg = _hierarchies[i].contractedPartitionedHypergraph();
      representative_hg.doParallelForAllNodes([&](const HypernodeID hn) {
        const HypernodeID coarse_hn = _hierarchies[i].mapToContractedHypergraph(hn);
        const PartitionID block = contracted_hg.partID(coarse_hn);
        ASSERT(block != kInvalidPartition && block < representative_hg.k());
        representative_hg.setOnlyNodePart(hn, block);
      });
      representative_hg.initializePartition(_task_group_id);

      ASSERT(metrics::objective(representative_hg, _context.partition.objective) ==
             metrics::objective(contracted_hg, _context.partition.objective),
             V(metrics::objective(representative_hg, _context.partition.objective)) <<
             V(metrics::objective(contracted_hg, _context.partition.objective)));
      ASSERT(metrics::imbalance(representative_hg, _context) ==
             metrics::imbalance(contracted_hg, _context),
             V(metrics::imbalance(representative_hg, _context)) <<
             V(metrics::imbalance(contracted_hg, _context)));
      utils::Timer::instance().stop_timer("projecting_partition");

      LOG << "Projecting partition" << V(i) << V(representative_hg.initialNumNodes()) << V(metrics::km1(contracted_hg)) << V(metrics::km1(representative_hg)) << V(current_metrics.km1) << "now refine";


      // Refinement
      refine(representative_hg, label_propagation, fm_refiner, current_metrics);

      // Update Progress Bar
      uncontraction_progress.setObjective(current_metrics.getMetric(_context.partition.mode, _context.partition.objective));
      uncontraction_progress += representative_hg.initialNumNodes() - contracted_hg.initialNumNodes();
    }

    // If we reach the original hypergraph and partition is imbalanced, we try to rebalance it
    if ( _top_level && metrics::imbalance(_partitioned_hg, _context) > _context.partition.epsilon) {
      const HyperedgeWeight quality_before = current_metrics.getMetric(
        kahypar::Mode::direct_kway, _context.partition.objective);
      if ( _context.partition.verbose_output ) {
        LOG << RED << "Partition is imbalanced (Current Imbalance:"
            << metrics::imbalance(_partitioned_hg, _context) << ") ->"
            << "Rebalancer is activated" << END;
      }

      utils::Timer::instance().start_timer("rebalance", "Rebalance");
      if ( _context.partition.objective == kahypar::Objective::km1 ) {
        Km1Rebalancer rebalancer(_partitioned_hg, _context, _task_group_id);
        rebalancer.rebalance(current_metrics);
      } else if ( _context.partition.objective == kahypar::Objective::cut ) {
        CutRebalancer rebalancer(_partitioned_hg, _context, _task_group_id);
        rebalancer.rebalance(current_metrics);
      }
      utils::Timer::instance().stop_timer("rebalance");

      const HyperedgeWeight quality_after = current_metrics.getMetric(
        kahypar::Mode::direct_kway, _context.partition.objective);
      if ( _context.partition.verbose_output ) {
        const HyperedgeWeight quality_delta = quality_after - quality_before;
        if ( quality_delta > 0 ) {
          LOG << RED << "Rebalancer worsen solution quality by" << quality_delta
              << "(Current Imbalance:" << metrics::imbalance(_partitioned_hg, _context) << ")" << END;
        } else {
          LOG << GREEN << "Rebalancer improves solution quality by" << abs(quality_delta)
              << "(Current Imbalance:" << metrics::imbalance(_partitioned_hg, _context) << ")" << END;
        }
      }
    }

    ASSERT(metrics::objective(_partitioned_hg, _context.partition.objective) ==
           current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
           V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective)) <<
           V(metrics::objective(_partitioned_hg, _context.partition.objective)));
    return std::move(_partitioned_hg);
  }

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
              refinement::MultiTryKWayFM& fm_refiner,
              kahypar::Metrics& current_metrics) {

    if ( label_propagation ) {
      utils::Timer::instance().start_timer("initialize_lp_refiner", "Initialize LP IRefiner");
      label_propagation->initialize(partitioned_hypergraph);
      utils::Timer::instance().stop_timer("initialize_lp_refiner");

      utils::Timer::instance().start_timer("label_propagation", "Label Propagation");
      label_propagation->refine(partitioned_hypergraph, current_metrics);
      utils::Timer::instance().stop_timer("label_propagation");
    }

    fm_refiner.refine(partitioned_hypergraph, current_metrics);
  }

  bool _is_finalized;
  Hypergraph& _hg;
  PartitionedHypergraph _partitioned_hg;
  const Context& _context;
  const TaskGroupID _task_group_id;
  const bool _top_level;
  parallel::scalable_vector<Hierarchy> _hierarchies;
};
}  // namespace mt_kahypar

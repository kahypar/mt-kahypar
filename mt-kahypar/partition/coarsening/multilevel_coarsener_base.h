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
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/progress_bar.h"
#include "mt-kahypar/utils/stats.h"

#include "mt-kahypar/partition/refinement/fm/multitry_kway_fm.h"  // for now. include it in registries later

namespace mt_kahypar {
template <typename TypeTraits>
class MultilevelCoarsenerBase {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using PartitionedHyperGraph = typename TypeTraits::PartitionedHyperGraph;
  using TBB = typename TypeTraits::TBB;

  using Refiner = IRefinerT<TypeTraits>;

  static constexpr bool debug = true;

  class Hierarchy {

   public:
    explicit Hierarchy(HyperGraph&& contracted_hypergraph,
                       parallel::scalable_vector<HypernodeID>&& communities) :
      _representative_hypergraph(nullptr),
      _contracted_hypergraph(std::move(contracted_hypergraph)),
      _contracted_partitioned_hypergraph(),
      _communities(std::move(communities)) { }

    void setRepresentativeHypergraph(PartitionedHyperGraph* representative_hypergraph) {
      _representative_hypergraph = representative_hypergraph;
    }

    PartitionedHyperGraph& representativeHypergraph() {
      ASSERT(_representative_hypergraph);
      return *_representative_hypergraph;
    }

    HyperGraph& contractedHypergraph() {
      return _contracted_hypergraph;
    }

    PartitionedHyperGraph& contractedPartitionedHypergraph() {
      return _contracted_partitioned_hypergraph;
    }

    const HyperGraph& contractedHypergraph() const {
      return _contracted_hypergraph;
    }

    // ! Maps a global vertex id of the representative hypergraph
    // ! to its global vertex id in the contracted hypergraph
    HypernodeID mapToContractedHypergraph(const HypernodeID hn) const {
      ASSERT(_representative_hypergraph);
      const HypernodeID original_id = _representative_hypergraph->originalNodeID(hn);
      ASSERT(original_id < _communities.size());
      return _communities[original_id];
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
    PartitionedHyperGraph* _representative_hypergraph;
    // ! Contracted Hypergraph
    HyperGraph _contracted_hypergraph;
    // ! Partitioned Hypergraph
    PartitionedHyperGraph _contracted_partitioned_hypergraph;
    // ! Defines the communities that are contracted
    // ! in the coarse hypergraph
    parallel::scalable_vector<HypernodeID> _communities;
  };

 public:
  MultilevelCoarsenerBase(HyperGraph& hypergraph, const Context& context, const TaskGroupID task_group_id) :
    uncontraction_progress(hypergraph.initialNumNodes(),
                           std::numeric_limits<HyperedgeWeight>::max(),
                           context.partition.verbose_output && context.partition.enable_progress_bar),
    _is_finalized(false),
    _hg(hypergraph),
    _partitioned_hg(),
    _context(context),
    _task_group_id(task_group_id),
    _hierarchies()
  {
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

  HyperGraph& coarsestHypergraph() {
    if ( _hierarchies.empty() ) {
      return _hg;
    } else {
      return _hierarchies.back().contractedHypergraph();
    }
  }

  PartitionedHyperGraph& coarsestPartitionedHypergraph() {
    ASSERT(_is_finalized);
    if ( _hierarchies.empty() ) {
      return _partitioned_hg;
    } else {
      return _hierarchies.back().contractedPartitionedHypergraph();
    }
  }

  void finalize() {
    utils::Timer::instance().start_timer("finalize_multilevel_hierarchy", "Finalize Multilevel Hierarchy");
    // Construct partitioned hypergraphs parallel
    tbb::task_group group;
    // Construct top level partitioned hypergraph
    group.run([&] {
      _partitioned_hg = PartitionedHyperGraph(_context.partition.k, _task_group_id, _hg);
    });
    // Construct partitioned hypergraph for each coarsened hypergraph in the hierarchy
    for ( size_t i = 0; i < _hierarchies.size(); ++i ) {
      group.run([&, i] {
        _hierarchies[i].contractedPartitionedHypergraph() =
                PartitionedHyperGraph(_context.partition.k, _task_group_id, _hierarchies[i].contractedHypergraph());
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
    HyperGraph& current_hg = coarsestHypergraph();
    ASSERT(current_hg.initialNumNodes() == communities.size());
    HyperGraph contracted_hg = current_hg.contract(communities, _task_group_id);
    _hierarchies.emplace_back(std::move(contracted_hg), std::move(communities));
  }

  PartitionedHyperGraph&& doUncoarsen(std::unique_ptr<Refiner>& label_propagation) {
    PartitionedHyperGraph& coarsest_hg = coarsestPartitionedHypergraph();
    initialize(coarsest_hg);

    refinement::MultiTryKWayFM fm_refiner(_context, _task_group_id, _hg.initialNumNodes(), _hg.initialNumEdges());

    // Refine Coarsest Partitioned Hypergraph
    refine(coarsest_hg, label_propagation, fm_refiner, current_metrics);

    for ( int i = _hierarchies.size() - 1; i >= 0; --i ) {
      // Project partition to next level finer hypergraph
      utils::Timer::instance().start_timer("projecting_partition", "Projecting Partition");
      PartitionedHyperGraph& representative_hg = _hierarchies[i].representativeHypergraph();
      PartitionedHyperGraph& contracted_hg = _hierarchies[i].contractedPartitionedHypergraph();

      representative_hg.doParallelForAllNodes(_task_group_id, [&](const HypernodeID hn) {
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


      // Refinement
      refine(representative_hg, label_propagation, fm_refiner, current_metrics);

      // Update Progress Bar
      uncontraction_progress.setObjective(current_metrics.getMetric(_context.partition.mode, _context.partition.objective));
      uncontraction_progress += representative_hg.initialNumNodes() - contracted_hg.initialNumNodes();
    }

    ASSERT(metrics::objective(_partitioned_hg, _context.partition.objective) ==
           current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
           V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective)) <<
           V(metrics::objective(_partitioned_hg, _context.partition.objective)));
    return std::move(_partitioned_hg);
  }

 protected:

  void computeMetrics(PartitionedHypergraph& phg) {
    tbb::parallel_invoke([&] { current_metrics.cut = metrics::hyperedgeCut(phg); }, [&] { current_metrics.km1 = metrics::km1(phg); });
    current_metrics.imbalance = metrics::imbalance(phg, _context);
  }

  void initialize(PartitionedHypergraph& current_hg) {
    computeMetrics(current_hg);
    int64_t num_nodes = current_hg.initialNumNodes();
    int64_t num_edges = current_hg.initialNumEdges();
    utils::Stats::instance().add_stat("initial_num_nodes", num_nodes);
    utils::Stats::instance().add_stat("initial_num_edges", num_edges);
    utils::Stats::instance().add_stat("initial_cut", current_metrics.cut);
    utils::Stats::instance().add_stat("initial_km1", current_metrics.km1);
    utils::Stats::instance().add_stat("initial_imbalance", current_metrics.imbalance);

    uncontraction_progress.setObjective(current_metrics.getMetric(_context.partition.mode, _context.partition.objective));
    uncontraction_progress += num_nodes;
  }

  kahypar::Metrics current_metrics;
  utils::ProgressBar uncontraction_progress;

  void refine(PartitionedHyperGraph& partitioned_hypergraph,
              std::unique_ptr<Refiner>& label_propagation,
              refinement::MultiTryKWayFM& fm_refiner,
              kahypar::Metrics& current_metrics) {
    if ( label_propagation ) {
      utils::Timer::instance().start_timer("initialize_lp_refiner", "Initialize LP Refiner");
      label_propagation->initialize(partitioned_hypergraph);
      utils::Timer::instance().stop_timer("initialize_lp_refiner");

      utils::Timer::instance().start_timer("label_propagation", "Label Propagation");
      label_propagation->refine(partitioned_hypergraph, {}, current_metrics);
      utils::Timer::instance().stop_timer("label_propagation");
    }


    if (_context.type == kahypar::ContextType::main) {
      fm_refiner.refine(partitioned_hypergraph, current_metrics);
    }

  }

  bool _is_finalized;
  HyperGraph& _hg;
  PartitionedHyperGraph _partitioned_hg;
  const Context& _context;
  const TaskGroupID _task_group_id;
  parallel::scalable_vector<Hierarchy> _hierarchies;
};
}  // namespace mt_kahypar

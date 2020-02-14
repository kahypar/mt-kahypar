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

namespace mt_kahypar {
template <typename TypeTraits>
class MultilevelCoarsenerBase {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using PartitionedHyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;
  using TBB = typename TypeTraits::TBB;

  using Refiner = IRefinerT<TypeTraits>;

  static constexpr bool debug = false;

  class Hierarchy {

   public:
    explicit Hierarchy(HyperGraph&& contracted_hypergraph,
                       const PartitionID k,
                       tbb::task_group& group,
                       parallel::scalable_vector<HypernodeID>&& communities,
                       parallel::scalable_vector<HypernodeID>&& mapping) :
      _representative_hypergraph(nullptr),
      _contracted_hypergraph(std::move(contracted_hypergraph)),
      _contracted_partitioned_hypergraph(),
      _communities(std::move(communities)),
      _mapping(std::move(mapping)) {
      ASSERT(_communities.size() == _mapping.size());
      // Construct partitioned hypergraph async
      group.run([&] {
        _contracted_partitioned_hypergraph = PartitionedHyperGraph(k, _contracted_hypergraph);
      });
    }

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
      return _mapping[_communities[original_id]];
    }

   private:
    // ! Hypergraph on the next finer level
    PartitionedHyperGraph* _representative_hypergraph;
    // ! Contracted Hypergraph
    HyperGraph _contracted_hypergraph;
    // ! Partitioned Hypergraph
    PartitionedHyperGraph _contracted_partitioned_hypergraph;
    // ! Defines the communities that are contracted
    // ! in the contracted hypergraph
    const parallel::scalable_vector<HypernodeID> _communities;
    // ! Mapping from community to original vertex id
    // ! in the contracted hypergraph
    const parallel::scalable_vector<HypernodeID> _mapping;
  };

 public:
  MultilevelCoarsenerBase(HyperGraph& hypergraph, const Context& context, const TaskGroupID task_group_id) :
    _hg(hypergraph),
    _partitioned_hg(),
    _context(context),
    _group(),
    _task_group_id(task_group_id),
    _hierarchies() {
    _group.run([&] {
      _partitioned_hg = PartitionedHyperGraph(_context.partition.k, _hg);
    });
  }

  MultilevelCoarsenerBase(const MultilevelCoarsenerBase&) = delete;
  MultilevelCoarsenerBase(MultilevelCoarsenerBase&&) = delete;
  MultilevelCoarsenerBase & operator= (const MultilevelCoarsenerBase &) = delete;
  MultilevelCoarsenerBase & operator= (MultilevelCoarsenerBase &&) = delete;

  virtual ~MultilevelCoarsenerBase() throw () {
    _group.wait();
  }

 protected:

  HypernodeID currentNumNodes() const {
    if ( _hierarchies.empty() ) {
      return _hg.initialNumNodes();
    } else {
      return _hierarchies.back().contractedHypergraph().initialNumNodes();
    }
  }

  HyperGraph& currentHypergraph() {
    if ( _hierarchies.empty() ) {
      return _hg;
    } else {
      return _hierarchies.back().contractedHypergraph();
    }
  }

  PartitionedHyperGraph& currentPartitionedHypergraph() {
    _group.wait();
    if ( _hierarchies.empty() ) {
      return _partitioned_hg;
    } else {
      return _hierarchies.back().contractedPartitionedHypergraph();
    }
  }

  void performMultilevelContraction(parallel::scalable_vector<HypernodeID>&& communities) {
    HyperGraph& current_hg = currentHypergraph();
    ASSERT(current_hg.initialNumNodes() == communities.size());
    auto contracted_hg = current_hg.contract(communities, _task_group_id);
    _hierarchies.emplace_back(std::move(contracted_hg.first), _context.partition.k, _group,
      std::move(communities), std::move(contracted_hg.second));
  }

  PartitionedHyperGraph&& doUncoarsen(std::unique_ptr<Refiner>& label_propagation) {
    const PartitionedHyperGraph& current_hg = currentPartitionedHypergraph();
    int64_t num_nodes = current_hg.initialNumNodes();
    int64_t num_edges = current_hg.initialNumEdges();
    HyperedgeWeight cut = 0;
    HyperedgeWeight km1 = 0;
    tbb::parallel_invoke([&] {
        // Cut metric
        cut = metrics::hyperedgeCut(current_hg);
      }, [&] {
        // Km1 metric
        km1 = metrics::km1(current_hg);
      });

    kahypar::Metrics current_metrics = { cut, km1, metrics::imbalance(current_hg, _context) };
    utils::Stats::instance().add_stat("initial_num_nodes", num_nodes);
    utils::Stats::instance().add_stat("initial_num_edges", num_edges);
    utils::Stats::instance().add_stat("initial_cut", current_metrics.cut);
    utils::Stats::instance().add_stat("initial_km1", current_metrics.km1);
    utils::Stats::instance().add_stat("initial_imbalance", current_metrics.imbalance);


    // We set the representative hypergraph of each hypergraph in the hierarchy
    // here, because due to resizing of the vector capacity it might become invalid
    // during creating the hierarchies
    _group.wait();
    if ( _hierarchies.size() > 0 ) {
      _hierarchies[0].setRepresentativeHypergraph(&_partitioned_hg);
      for ( size_t i = 1; i < _hierarchies.size(); ++i ) {
        _hierarchies[i].setRepresentativeHypergraph(&_hierarchies[i - 1].contractedPartitionedHypergraph());
      }
    }

    utils::ProgressBar uncontraction_progress(_hg.initialNumNodes(),
      _context.partition.objective == kahypar::Objective::km1 ? current_metrics.km1 : current_metrics.cut,
      _context.partition.verbose_output && _context.partition.enable_progress_bar);
    uncontraction_progress += num_nodes;

    for ( int i = _hierarchies.size() - 1; i >= 0; --i ) {
      // Project partition to next level finer hypergraph
      utils::Timer::instance().start_timer("projecting_partition", "Projecting Partition");
      PartitionedHyperGraph& representative_hg = _hierarchies[i].representativeHypergraph();
      PartitionedHyperGraph& contracted_hg = _hierarchies[i].contractedPartitionedHypergraph();
      tbb::parallel_for(0UL, representative_hg.initialNumNodes(), [&](const HypernodeID id) {
        const HypernodeID hn = representative_hg.globalNodeID(id);
        if ( representative_hg.nodeIsEnabled(hn) ) {
          const HypernodeID coarse_hn = _hierarchies[i].mapToContractedHypergraph(hn);
          const PartitionID block = contracted_hg.partID(coarse_hn);
          ASSERT(block != -1 && block < representative_hg.k());
          representative_hg.setNodePart(hn, block);
        }
      });
      representative_hg.initializeNumCutHyperedges(_task_group_id);
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
      utils::Timer::instance().start_timer("initialize_refiner", "Initialize Refiner");
      if ( label_propagation ) {
        label_propagation->initialize(representative_hg);
      }
      utils::Timer::instance().stop_timer("initialize_refiner");

      if ( label_propagation ) {
        label_propagation->refine(representative_hg, {}, current_metrics);
      }

      // Update Progress Bar
      uncontraction_progress.setObjective(
        _context.partition.objective == kahypar::Objective::km1 ?
        current_metrics.km1 : current_metrics.cut);
      uncontraction_progress += representative_hg.initialNumNodes() - contracted_hg.initialNumNodes();
    }

    ASSERT(metrics::objective(_partitioned_hg, _context.partition.objective) ==
           current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
           V(current_metrics.getMetric(kahypar::Mode::direct_kway, _context.partition.objective)) <<
           V(metrics::objective(_partitioned_hg, _context.partition.objective)));
    return std::move(_partitioned_hg);
  }

 protected:
  HyperGraph& _hg;
  PartitionedHyperGraph _partitioned_hg;
  const Context& _context;
  tbb::task_group _group;
  const TaskGroupID _task_group_id;
  parallel::scalable_vector<Hierarchy> _hierarchies;
};
}  // namespace mt_kahypar

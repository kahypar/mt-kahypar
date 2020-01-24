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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
template <typename TypeTraits>
class MultilevelCoarsenerBase {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
  using TBB = typename TypeTraits::TBB;

  using Memento = typename StreamingHyperGraph::Memento;
  using HypergraphPruner = HypergraphPrunerT<TypeTraits>;

  static constexpr bool debug = false;

  class Hierarchy {

   public:
    explicit Hierarchy(HyperGraph&& contracted_hypergraph,
                       parallel::scalable_vector<HypernodeID>&& communities,
                       parallel::scalable_vector<HypernodeID>&& mapping) :
      _representative_hypergraph(nullptr),
      _contracted_hypergraph(std::move(contracted_hypergraph)),
      _communities(std::move(communities)),
      _mapping(std::move(mapping)) {
      ASSERT(_communities.size() == _mapping.size());
    }

    void setRepresentativeHypergraph(HyperGraph* representative_hypergraph) {
      _representative_hypergraph = representative_hypergraph;
    }

    HyperGraph& representativeHypergraph() {
      ASSERT(_representative_hypergraph);
      return *_representative_hypergraph;
    }

    HyperGraph& contractedHypergraph() {
      return _contracted_hypergraph;
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
      const HypernodeID original_contracted_id = _mapping[_communities[original_id]];
      ASSERT(original_contracted_id < _contracted_hypergraph.initialNumNodes());
      return _contracted_hypergraph.globalNodeID(original_contracted_id);
    }

   private:
    HyperGraph* _representative_hypergraph;
    HyperGraph _contracted_hypergraph;
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
    _context(context),
    _task_group_id(task_group_id),
    _hierarchies() { }

  MultilevelCoarsenerBase(const MultilevelCoarsenerBase&) = delete;
  MultilevelCoarsenerBase(MultilevelCoarsenerBase&&) = delete;
  MultilevelCoarsenerBase & operator= (const MultilevelCoarsenerBase &) = delete;
  MultilevelCoarsenerBase & operator= (MultilevelCoarsenerBase &&) = delete;

  virtual ~MultilevelCoarsenerBase() throw () { }

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

  void performMultilevelContraction(parallel::scalable_vector<HypernodeID>&& communities) {
    HyperGraph& current_hg = currentHypergraph();
    ASSERT(current_hg.initialNumNodes() == communities.size());
    auto contracted_hg = current_hg.contract(communities, _task_group_id);
    _hierarchies.emplace_back(std::move(contracted_hg.first),
      std::move(communities), std::move(contracted_hg.second));
  }

  bool doUncoarsen(std::unique_ptr<IRefiner>&) {
    // We set the representative hypergraph of each hypergraph in the hierarchy
    // here, because due to resizing of the vector capacity it might become invalid
    // during creating the hierarchies
    if ( _hierarchies.size() > 0 ) {
      _hierarchies[0].setRepresentativeHypergraph(&_hg);
      for ( size_t i = 1; i < _hierarchies.size(); ++i ) {
        _hierarchies[i].setRepresentativeHypergraph(&_hierarchies[i - 1].contractedHypergraph());
      }
    }


    while ( !_hierarchies.empty() ) {
      // Project partition to next level finer hypergraph
      HyperGraph& representative_hg = _hierarchies.back().representativeHypergraph();
      HyperGraph& contracted_hg = _hierarchies.back().contractedHypergraph();
      tbb::parallel_for(0UL, representative_hg.initialNumNodes(), [&](const HypernodeID id) {
        const HypernodeID hn = representative_hg.globalNodeID(id);
        if ( representative_hg.nodeIsEnabled(hn) ) {
          const HypernodeID coarse_hn = _hierarchies.back().mapToContractedHypergraph(hn);
          const PartitionID block = contracted_hg.partID(coarse_hn);
          ASSERT(block != -1 && block < representative_hg.k());
          representative_hg.setNodePart(hn, block);
        }
      });
      representative_hg.updateGlobalPartInfos();
      representative_hg.initializeNumCutHyperedges();
      ASSERT(metrics::objective(representative_hg, _context.partition.objective) ==
             metrics::objective(contracted_hg, _context.partition.objective),
             V(metrics::objective(representative_hg, _context.partition.objective)) <<
             V(metrics::objective(contracted_hg, _context.partition.objective)));
      ASSERT(metrics::imbalance(representative_hg, _context) ==
             metrics::imbalance(contracted_hg, _context),
             V(metrics::imbalance(representative_hg, _context)) <<
             V(metrics::imbalance(contracted_hg, _context)));
      _hierarchies.pop_back();

      // TODO: Do some refinement stuff here
    }

    return true;
  }

 protected:
  HyperGraph& _hg;
  const Context& _context;
  const TaskGroupID _task_group_id;
  parallel::scalable_vector<Hierarchy> _hierarchies;
};
}  // namespace mt_kahypar

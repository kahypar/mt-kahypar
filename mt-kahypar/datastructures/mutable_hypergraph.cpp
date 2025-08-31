#include "mutable_hypergraph.h"

#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/datastructures/concurrent_bucket_map.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/memory_tree.h"

#include <tbb/parallel_reduce.h>
#include <tbb/parallel_sort.h>

namespace mt_kahypar::ds {


    /*!
    * This struct is used during multilevel coarsening to efficiently
    * detect parallel hyperedges.
    */
    struct ContractedHyperedgeInformation {
        HyperedgeID he = kInvalidHyperedge;
        size_t hash = kEdgeHashSeed;
        size_t size = std::numeric_limits<size_t>::max();
        bool valid = false;
    };

    /*!
     * Contracts a given community structure. All vertices with the same label
     * are collapsed into the same vertex. The resulting single-pin and parallel
     * hyperedges are removed from the contracted graph. The function returns
     * the contracted hypergraph and a mapping which specifies a mapping from
     * community label (given in 'communities') to a vertex in the coarse hypergraph.
     *
     * \param communities Community structure that should be contracted
     */
    MutableHypergraph MutableHypergraph::contract_transform(parallel::scalable_vector<HypernodeID>& communities, bool deterministic) {

      ASSERT(communities.size() >= _num_hypernodes);

      std::vector<HypernodeID> static_to_mut_hn;
      std::vector<HyperedgeID> static_to_mut_he;
      std::vector<HypernodeID> deleted_hn;
      std::vector<HyperedgeID> deleted_he;

      StaticHypergraph static_hg = toStaticHypergraph(&static_to_mut_hn, &static_to_mut_he, &deleted_hn, &deleted_he);
      parallel::scalable_vector<HypernodeID> static_communities = toStaticCommunity(&communities);

      StaticHypergraph contracted_static_hg = static_hg.contract(static_communities, deterministic);

      MutableHypergraph contracted_hg = fromStaticHypergraph(contracted_static_hg, &static_to_mut_hn, &static_to_mut_he, &deleted_hn, &deleted_he);
      updateFromStaticCommunity(static_communities, &communities, &static_to_mut_hn);

      return contracted_hg;
    }

    /*!
     * Contracts a given community structure. All vertices with the same label
     * are collapsed into the same vertex. The resulting single-pin and parallel
     * hyperedges are removed from the contracted graph. The function returns
     * the contracted hypergraph and a mapping which specifies a mapping from
     * community label (given in 'communities') to a vertex in the coarse hypergraph.
     *
     * \param communities Community structure that should be contracted
     */
    MutableHypergraph MutableHypergraph::contract(parallel::scalable_vector<HypernodeID>& communities, bool deterministic) {
      (void) deterministic;

      ASSERT(communities.size() >= _num_hypernodes);


      // filter all invalid communities
      std::vector<HypernodeID> unique_communities;
      for (size_t i = 0; i < communities.size(); ++i) {
        if (communities[i] != kInvalidHypernode && !_hypernodes[i].is_deleted() && !_hypernodes[i].isDisabled()) {
          unique_communities.push_back(communities[i]);
        }
      }

      std::sort(unique_communities.begin(), unique_communities.end());
      auto last = std::unique(unique_communities.begin(), unique_communities.end());
      unique_communities.erase(last, unique_communities.end());

      std::vector<HypernodeID> community_map(unique_communities.back() + 1);
      for (size_t i = 0; i < unique_communities.size(); ++i) {
        if (unique_communities[i] == kInvalidHypernode || (_hypernodes[unique_communities[i]].is_deleted() || _hypernodes[unique_communities[i]].isDisabled())) {
          continue;
        }
        community_map[unique_communities[i]] = i;
      }

      //remap communities [5, 2, 5] -> [0, 1, 0]
      std::unordered_map<size_t, size_t> mapping; // value -> remapped index
      for (size_t i = 0; i < communities.size(); ++i) {
        if (communities[i] == kInvalidHypernode || _hypernodes[i].is_deleted() || _hypernodes[i].isDisabled()) {
          continue;
        }
        communities[i] = community_map[communities[i]];
      }

      MutableHypergraph hypergraph;

      hypergraph._num_hypernodes = unique_communities.size();
      hypergraph._hypernodes.resize(hypergraph._num_hypernodes + 1); // +1 for sentinel
      for (size_t hn = 0; hn < communities.size(); ++hn) {
        if (communities[hn] == kInvalidHypernode || _hypernodes[hn].is_deleted() || _hypernodes[hn].isDisabled()) {
          continue;
        }
        HypernodeID coarse_hn = communities[hn];
        if (!hypergraph.nodeIsEnabled(coarse_hn)) {
          hypergraph.hypernode(coarse_hn).enable();
          hypergraph.hypernode(coarse_hn).setWeight(0);
        }
        hypergraph.hypernode(coarse_hn).setWeight(hypergraph.hypernode(coarse_hn).weight() + _hypernodes[hn].weight());
        for (const HyperedgeID he : _hypernodes[hn].incidentEdges()) {
          if (std::find(hypergraph.incidentEdges(coarse_hn).begin(), hypergraph.incidentEdges(coarse_hn).end(), he) != hypergraph.incidentEdges(coarse_hn).end()) {
            continue;
          }
          hypergraph._hypernodes[coarse_hn].addNet(he);
        }
      }

      std::vector<std::vector<HypernodeID>> new_pins;

      for (HypernodeID hn : hypergraph.nodes()) {
        for (const HyperedgeID he : hypergraph.incidentEdges(hn)) {
          if (_hyperedges[he].is_deleted() || _hyperedges[he].isDisabled()) {
            continue;
          }
          if (he >= new_pins.size()) {
            new_pins.resize(he + 1);
          }
          new_pins[he].push_back(hn);
        }
      }

      std::vector<std::vector<HypernodeID>> filtered_pins;
      std::vector<HyperedgeWeight> filtered_weights;


      //sort pins and remove duplicates
      for (size_t he = 0; he < new_pins.size(); ++he) {
        std::sort(new_pins[he].begin(), new_pins[he].end());
        auto last = std::unique(new_pins[he].begin(), new_pins[he].end());
        new_pins[he].erase(last, new_pins[he].end());
        if (new_pins[he].size() <= 1) {
          continue;
        } else {
          auto it = std::find_if(filtered_pins.begin(), filtered_pins.end(),
                             [&](const std::vector<HypernodeID>& pins) {
                               return pins == new_pins[he];
                             });
          if (it == filtered_pins.end()) {
            filtered_pins.push_back(new_pins[he]);
            filtered_weights.push_back(1);
          } else {
            //parallel hyperedge -> increase weight of existing hyperedge
            size_t index = std::distance(filtered_pins.begin(), it);
            filtered_weights[index] += 1;
          }
        }
      }

      //delete all incident nets
      for (HypernodeID hn : hypergraph.nodes()) {
        for (const HyperedgeID he : hypergraph.incidentEdges(hn)) {
          hypergraph.hypernode(hn).deleteNet(he);
        }
      }

      hypergraph._hyperedges.resize(filtered_pins.size());
      hypergraph._num_hyperedges = hypergraph._hyperedges.size();
      for (size_t he = 0; he < filtered_pins.size(); ++he) {
        hypergraph._hyperedges[he] = Hyperedge();
        hypergraph._hyperedges[he].enable();
        hypergraph._hyperedges[he].setWeight(filtered_weights[he]);
        for (const HypernodeID hn : filtered_pins[he]) {
          hypergraph.addPin(he, hn);
        }
      }

      //TODO community ids
//      auto assign_communities = [&] {
//          hypergraph._community_ids.resize(num_hypernodes, 0);
//          doParallelForAllNodes([&](HypernodeID fine_hn) {
//              hypergraph.setCommunityID(map_to_coarse_hypergraph(fine_hn), communityID(fine_hn));
//          });
//      };
      hypergraph._community_ids.resize(hypergraph._num_hypernodes, 0);
      for (size_t hn : nodes()) {
        hypergraph.setCommunityID(communities[hn], _community_ids[hn]);
      }

      hypergraph._total_weight = totalWeight();
      return hypergraph;
    }


    parallel::scalable_vector<HypernodeID> MutableHypergraph::toStaticCommunity(parallel::scalable_vector<HypernodeID>* mutable_communities) {
      parallel::scalable_vector<HypernodeID> static_communities = parallel::scalable_vector<HypernodeID>(_num_hypernodes, kInvalidHypernode);

      // create mapping from mutable hypernode ids to static hypernode ids
      HypernodeID static_hn_id = 0;
      for (size_t hn = 0; hn < mutable_communities->size(); ++hn) {
        if (_hypernodes[hn].is_deleted()) {
          continue;
        }
        static_communities[static_hn_id] = mutable_communities->at(hn);
        ++static_hn_id;
      }
      return static_communities;
    }

    // ! Convert from static community
    void MutableHypergraph::updateFromStaticCommunity(const parallel::scalable_vector<HypernodeID>& static_communities, parallel::scalable_vector<HypernodeID>* mutable_communities,
                                                               std::vector<HypernodeID>* static_to_mut_hn) {
      for (size_t static_hn = 0; static_hn < static_communities.size(); ++static_hn) {
        HypernodeID mut_hn = static_to_mut_hn->at(static_hn);
        mutable_communities->at(mut_hn) = static_communities[static_hn];
      }
    }


    // ! Copy static hypergraph in parallel
    MutableHypergraph MutableHypergraph::copy(parallel_tag_t) const {
      // parallel copy is not supported in mutable hypergraph
      return copy();
    }

    // ! Copy static hypergraph sequential
    MutableHypergraph MutableHypergraph::copy() const {
      MutableHypergraph hypergraph;

      hypergraph._num_hypernodes = _num_hypernodes;
      hypergraph._num_removed_hypernodes = _num_removed_hypernodes;
      hypergraph._num_hyperedges = _num_hyperedges;
      hypergraph._num_removed_hyperedges = _num_removed_hyperedges;
      hypergraph._max_edge_size = _max_edge_size;
      hypergraph._num_pins = _num_pins;
      hypergraph._total_degree = _total_degree;
      hypergraph._total_weight = _total_weight;

      hypergraph._hypernodes.resize(_hypernodes.size());
      for (HypernodeID hn : nodes()) {
        hypergraph._hypernodes[hn] = Hypernode(true);
        hypergraph._hypernodes[hn].setWeight(_hypernodes[hn].weight());
        for (const HyperedgeID he : _hypernodes[hn].incidentEdges()) {
          hypergraph._hypernodes[hn].addNet(he);
        }
        if (_hypernodes[hn].isDisabled()) {
          hypergraph._hypernodes[hn].disable();
        }
        if (_hypernodes[hn].is_deleted()){
          hypergraph._hypernodes[hn].mark_deleted();
        }
      }

      hypergraph._hyperedges.resize(_hyperedges.size());
      for (HyperedgeID he : edges()) {
        hypergraph._hyperedges[he] = Hyperedge();
        if (!_hyperedges[he].isDisabled()) {
          hypergraph._hyperedges[he].enable();
        }
        hypergraph._hyperedges[he].setWeight(_hyperedges[he].weight());
        for (const HypernodeID hn : _hyperedges[he].pins()) {
          hypergraph._hyperedges[he].addPin(hn);
        }
        if (_hyperedges[he].is_deleted()){
          hypergraph._hyperedges[he].mark_deleted();
        }
      }

      hypergraph._community_ids = _community_ids;
      hypergraph.addFixedVertexSupport(_fixed_vertices.copy());

      return hypergraph;
    }

    StaticHypergraph MutableHypergraph::toStaticHypergraph(std::vector<HypernodeID>* static_to_mut_hn,
                                                           std::vector<HyperedgeID>* static_to_mut_he, std::vector<HypernodeID>* deleted_hn, std::vector<HyperedgeID>* deleted_he) {
      StaticHypergraph hypergraph;

      hypergraph._num_hypernodes = _num_hypernodes;
      hypergraph._num_removed_hypernodes = _num_removed_hypernodes;
      hypergraph._num_hyperedges = _num_hyperedges;
      hypergraph._num_removed_hyperedges = _num_removed_hyperedges;
      hypergraph._max_edge_size = _max_edge_size;
      hypergraph._num_pins = _num_pins;
      hypergraph._total_degree = _total_degree;
      hypergraph._total_weight = _total_weight;

      hypergraph._community_ids = _community_ids;
      //TODO add fixed vertex support

      // create incidence array and incident nets
      hypergraph._incidence_array.resize(_num_pins);
      hypergraph._hyperedges.resize(_num_hyperedges + 1);
      size_t incidence_array_pos = 0;

      // resize mapping vectors
      static_to_mut_he->resize(_hyperedges.size());
      size_t mutable_he_index = -1;
      size_t static_he_index = -1;

      for (size_t he = 0; he < _num_hyperedges; ++he) {
        ++mutable_he_index;
        if (_hyperedges[he].is_deleted()) {
          deleted_he->push_back(mutable_he_index);
          continue;
        }
        ++static_he_index;
        static_to_mut_he->at(static_he_index) = mutable_he_index;
        auto hyperedge_test = StaticHypergraph::Hyperedge();
        hypergraph._hyperedges[he] = hyperedge_test;
        hypergraph._hyperedges[he].enable();
        // temporarly enable mutable hypernode to allow copy
        bool is_disabled = !edgeIsEnabled(he);
        if (is_disabled) _hyperedges[he].enable();
        hypergraph._hyperedges[he].setWeight(_hyperedges[he].weight());
        hypergraph._hyperedges[he].setSize(_hyperedges[he].size());
        hypergraph._hyperedges[he].setFirstEntry(incidence_array_pos);

        for (const HypernodeID hn : _hyperedges[he].pins()) {
          hypergraph._incidence_array[incidence_array_pos++] = hn;
        }
        if (is_disabled) {
          hypergraph._hyperedges[he].disable();
          _hyperedges[he].disable();
        }
      }

      hypergraph._incident_nets.resize(_total_degree);
      hypergraph._hypernodes.resize(_num_hypernodes + 1);
      size_t incident_nets_pos = 0;

      static_to_mut_hn->resize(_hypernodes.size());
      size_t mutable_hn_index = -1;
      size_t static_hn_index = -1;

      for (size_t hn = 0; hn < _num_hypernodes; ++hn) {
        ++mutable_hn_index;
        if (_hypernodes[hn].is_deleted()) {
          deleted_hn->push_back(mutable_hn_index);
          continue;
        }
        ++static_hn_index;
        static_to_mut_hn->at(static_hn_index) = mutable_hn_index;
        hypergraph._hypernodes[hn] = StaticHypergraph::Hypernode();
        hypergraph._hypernodes[hn].enable();
        // temporarly enable mutable hypernode to allow copy
        bool is_disabled = !nodeIsEnabled(hn);
        if (is_disabled) _hypernodes[hn].enable();
        hypergraph._hypernodes[hn].setWeight(_hypernodes[hn].weight());
        hypergraph._hypernodes[hn].setSize(_hypernodes[hn].size());
        hypergraph._hypernodes[hn].setFirstEntry(incident_nets_pos);
        for (const HyperedgeID he : _hypernodes[hn].incidentEdges()) {
          hypergraph._incident_nets[incident_nets_pos++] = he;
        }
        if (is_disabled) {
          hypergraph._hypernodes[hn].disable();
          _hypernodes[hn].disable();
        }
      }

      // TODO Add sentinels ?
      hypergraph._hypernodes.back() = StaticHypergraph::Hypernode(hypergraph._incident_nets.size());
      hypergraph._hyperedges.back() = StaticHypergraph::Hyperedge(hypergraph._incidence_array.size());

      return hypergraph;
    }

    MutableHypergraph MutableHypergraph::fromStaticHypergraph(const StaticHypergraph &static_hg, std::vector<HypernodeID>* static_to_mut_hn,
                                                              std::vector<HyperedgeID>* static_to_mut_he, std::vector<HypernodeID>* deleted_hn, std::vector<HyperedgeID>* deleted_he) {
      MutableHypergraph hypergraph;

      hypergraph._num_hypernodes = static_hg._num_hypernodes;
      hypergraph._num_removed_hypernodes = static_hg._num_removed_hypernodes;
      hypergraph._num_hyperedges = static_hg._num_hyperedges;
      hypergraph._num_removed_hyperedges = static_hg._num_removed_hyperedges;
      hypergraph._max_edge_size = static_hg._max_edge_size;
      hypergraph._num_pins = static_hg._num_pins;
      hypergraph._total_degree = static_hg._total_degree;
      hypergraph._total_weight = static_hg._total_weight;

      hypergraph._hypernodes.resize(static_hg._num_hypernodes);

      for (HypernodeID hn : static_hg.nodes()) {
        const StaticHypergraph::Hypernode& node = static_hg.hypernode(hn);
        const size_t mutable_hn_index = static_to_mut_hn->at(hn);
        hypergraph._hypernodes[mutable_hn_index] = MutableHypergraph::Hypernode(static_hg.nodeIsEnabled(hn));
        for (HyperedgeID he : static_hg.incidentEdges(hn)) {
          hypergraph._hypernodes[mutable_hn_index].addNet(static_to_mut_he->at(he));
        }
        hypergraph._hypernodes[mutable_hn_index].setWeight(node.weight());
      }

      for (HypernodeID hn : *deleted_hn) {
        hypergraph._hypernodes[hn] = MutableHypergraph::Hypernode();
        hypergraph._hypernodes[hn].mark_deleted();
      }

      hypergraph._hyperedges.resize(static_hg._num_hyperedges);
      for (HyperedgeID he : static_hg.edges()) {
        const StaticHypergraph::Hyperedge& edge = static_hg.hyperedge(he);
        const size_t mutable_he_index = static_to_mut_he->at(he);
        hypergraph._hyperedges[mutable_he_index] = MutableHypergraph::Hyperedge();
        hypergraph._hyperedges[mutable_he_index].enable();
        hypergraph._hyperedges[mutable_he_index].setWeight(edge.weight());
        for (HypernodeID hn : static_hg.pins(he)) {
          hypergraph._hyperedges[mutable_he_index].addPin(hn);
        }
        if (!static_hg.edgeIsEnabled(he)) {
          hypergraph._hyperedges[mutable_he_index].disable();
        }
      }

      for (HyperedgeID he : *deleted_he) {
        hypergraph._hyperedges[he] = MutableHypergraph::Hyperedge();
        hypergraph._hyperedges[he].mark_deleted();
      }

      hypergraph._community_ids = static_hg._community_ids;

      // Add sentinels
      hypergraph._hypernodes.push_back(MutableHypergraph::Hypernode(static_cast<size_t>(hypergraph._total_degree)));
      hypergraph._hyperedges.push_back(MutableHypergraph::Hyperedge(hypergraph._num_pins));

      return hypergraph;
    }

    void MutableHypergraph::memoryConsumption(utils::MemoryTreeNode* parent) const {
      ASSERT(parent);
      parent->addChild("Hypernodes", sizeof(Hypernode) * _hypernodes.size());
      parent->addChild("Hyperedges", sizeof(Hyperedge) * _hyperedges.size());
      parent->addChild("Communities", sizeof(PartitionID) * _community_ids.capacity());
      if ( hasFixedVertices() ) {
        parent->addChild("Fixed Vertex Support", _fixed_vertices.size_in_bytes());
      }
    }

    // ! Computes the total node weight of the hypergraph
    void MutableHypergraph::computeAndSetTotalNodeWeight(parallel_tag_t) {
      _total_weight = tbb::parallel_reduce(tbb::blocked_range<HypernodeID>(ID(0), _num_hypernodes), 0,
                                           [this](const tbb::blocked_range<HypernodeID>& range, HypernodeWeight init) {
                                               HypernodeWeight weight = init;
                                               for (HypernodeID hn = range.begin(); hn < range.end(); ++hn) {
                                                 if (nodeIsEnabled(hn)) {
                                                   weight += this->_hypernodes[hn].weight();
                                                 }
                                               }
                                               return weight;
                                           }, std::plus<>());
    }

} // namespace

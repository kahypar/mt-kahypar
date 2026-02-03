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
    MutableHypergraph MutableHypergraph::contract(parallel::scalable_vector<HypernodeID>& communities, bool deterministic) {

      ASSERT(communities.size() == _num_hypernodes);

      StaticHypergraph static_hg = toStaticHypergraph();

      StaticHypergraph contracted_static_hg = static_hg.contract(communities, deterministic);

      MutableHypergraph contracted_hg = fromStaticHypergraph(contracted_static_hg);

      return contracted_hg;
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

    StaticHypergraph MutableHypergraph::toStaticHypergraph() {
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
      hypergraph._hyperedges.resize(_num_hyperedges);
      size_t incidence_array_pos = 0;

      for (HyperedgeID he : edges()) {
        auto hyperedge_test = StaticHypergraph::Hyperedge();
        hypergraph._hyperedges[he] = hyperedge_test;
        if (!hyperedge(he).isDisabled()) {
          hypergraph._hyperedges[he].enable();
        }
        hypergraph._hyperedges[he].setWeight(_hyperedges[he].weight());
        hypergraph._hyperedges[he].setSize(_hyperedges[he].size());
        hypergraph._hyperedges[he].setFirstEntry(incidence_array_pos);

        for (const HypernodeID hn : _hyperedges[he].pins()) {
          hypergraph._incidence_array[incidence_array_pos++] = hn;
        }
      }

      hypergraph._incident_nets.resize(_total_degree);
      hypergraph._hypernodes.resize(_num_hypernodes);
      size_t incident_nets_pos = 0;

      for (HypernodeID hn : nodes()) {
        hypergraph._hypernodes[hn] = StaticHypergraph::Hypernode();
        if (nodeIsEnabled(hn)) {
          hypergraph._hypernodes[hn].enable();
        }
        hypergraph._hypernodes[hn].setWeight(_hypernodes[hn].weight());
        hypergraph._hypernodes[hn].setSize(_hypernodes[hn].size());
        hypergraph._hypernodes[hn].setFirstEntry(incident_nets_pos);
        for (const HyperedgeID he : _hypernodes[hn].incidentEdges()) {
          hypergraph._incident_nets[incident_nets_pos++] = he;
        }
      }

      return hypergraph;
    }

    MutableHypergraph MutableHypergraph::fromStaticHypergraph(const StaticHypergraph &static_hg) {
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
        hypergraph._hypernodes[hn] = MutableHypergraph::Hypernode(static_hg.nodeIsEnabled(hn));
        for (HyperedgeID he : static_hg.incidentEdges(hn)) {
          hypergraph._hypernodes[hn].addNet(he);
        }
        hypergraph._hypernodes[hn].setWeight(node.weight());
      }

      hypergraph._hyperedges.resize(static_hg._num_hyperedges);
      for (HyperedgeID he : static_hg.edges()) {
        const StaticHypergraph::Hyperedge& edge = static_hg.hyperedge(he);
        hypergraph._hyperedges[he] = MutableHypergraph::Hyperedge();
        hypergraph._hyperedges[he].enable();
        hypergraph._hyperedges[he].setWeight(edge.weight());
        for (HypernodeID hn : static_hg.pins(he)) {
          hypergraph._hyperedges[he].addPin(hn);
        }
        if (!static_hg.edgeIsEnabled(he)) {
          hypergraph._hyperedges[he].disable();
        }
      }
      hypergraph._community_ids = static_hg._community_ids;

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

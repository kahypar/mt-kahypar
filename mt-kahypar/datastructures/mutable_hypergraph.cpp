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

      MutableHypergraph contracted_hg = fromStaticHypergraphSimple(contracted_static_hg);
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

      if ( !_tmp_contraction_buffer ) {
        allocateTmpContractionBuffer();
      }

      // Auxiliary buffers - reused during multilevel hierarchy to prevent expensive allocations
      // std::vector<size_t>& mapping = _tmp_contraction_buffer->mapping;
      // std::vector<Hypernode>& tmp_hypernodes = _tmp_contraction_buffer->tmp_hypernodes;
      // std::vector<parallel::IntegralAtomicWrapper<size_t>>& tmp_num_incident_nets =
      //         _tmp_contraction_buffer->tmp_num_incident_nets;
      // std::vector<parallel::IntegralAtomicWrapper<HypernodeWeight>>& hn_weights =
      //         _tmp_contraction_buffer->hn_weights;
      std::vector<Hyperedge>& tmp_hyperedges = _tmp_contraction_buffer->tmp_hyperedges;
      // std::vector<size_t>& he_sizes = _tmp_contraction_buffer->he_sizes;
      std::vector<size_t>& valid_hyperedges = _tmp_contraction_buffer->valid_hyperedges;

      ASSERT(static_cast<size_t>(_num_hyperedges) <= tmp_hyperedges.size());
      ASSERT(static_cast<size_t>(_num_hyperedges) <= valid_hyperedges.size());

      // #################### STAGE 1 ####################
      // The communities are mapped to a contiguous range [0, k]
      // The order of communities is preserved, i.e. [5, 2, 5] -> [1, 0, 1]
      // Invalid nodes are mapped to kInvalidHypernode

      // filter all invalid communities
      std::vector<HypernodeID> unique_communities;
      for (size_t i = 0; i < communities.size(); ++i) {
        if (communities[i] != kInvalidHypernode && !_hypernodes[i].is_deleted() && !_hypernodes[i].isDisabled()) {
          unique_communities.push_back(communities[i]);
        }
      }

      //sort and remove duplicates
      std::sort(unique_communities.begin(), unique_communities.end());
      auto last = std::unique(unique_communities.begin(), unique_communities.end());
      unique_communities.erase(last, unique_communities.end());

      //create mapping from old community id to new community id
      std::vector<HypernodeID> community_map(unique_communities.back() + 1, kInvalidHypernode);
      for (size_t i = 0; i < unique_communities.size(); ++i) {
        community_map[unique_communities[i]] = i;
      }

      //remap communities
      for (size_t i = 0; i < communities.size(); ++i) {
        if (communities[i] == kInvalidHypernode || _hypernodes[i].is_deleted() || _hypernodes[i].isDisabled()) {
          continue;
        }
        communities[i] = community_map[communities[i]];
      }


      // #################### STAGE 2 ####################
      // Create a node for each community with summed weight and union of incident nets

      auto cs2 = [](const HypernodeID x) { return x * x; };
      ConcurrentBucketMap<ContractedHyperedgeInformation> hyperedge_hash_map;
      hyperedge_hash_map.reserve_for_estimated_number_of_insertions(_num_hyperedges);

      for (HyperedgeID he : edges()) {
        if (edgeIsEnabled(he)) {
          // Copy hyperedge and pins to temporary buffer
          const Hyperedge &e = _hyperedges[he];
          ASSERT(static_cast<size_t>(he) < tmp_hyperedges.size());
          tmp_hyperedges[he] = e;
          valid_hyperedges[he] = 1;

          // Map pins to vertex ids in coarse graph
          std::vector<HypernodeID> pins;
          for (HypernodeID hn : e.pins()) {
            if (communities[hn] == kInvalidHypernode || _hypernodes[hn].is_deleted() || _hypernodes[hn].isDisabled()) {
              continue;
            }
            pins.push_back(communities[hn]);
          }

          // Remove duplicates and disabled vertices
          std::sort(pins.begin(), pins.end());
          auto last = std::unique(pins.begin(), pins.end());
          pins.erase(last, pins.end());

          // Update (size of) hyperedge in temporary hyperedge buffer
          const size_t contracted_size = pins.size();
          tmp_hyperedges[he] = Hyperedge(pins, e.weight());


          if (contracted_size > 1) {
            // Compute hash of contracted hyperedge
            size_t footprint = kEdgeHashSeed;
            for (size_t pos = 0; pos < contracted_size; ++pos) {
              footprint += cs2(pins[pos]);
            }
            hyperedge_hash_map.insert(footprint,
                                      ContractedHyperedgeInformation{he, footprint, contracted_size, true});
          } else {
            // Hyperedge becomes a single-pin hyperedge
            valid_hyperedges[he] = 0;
            tmp_hyperedges[he].disable();
          }
        } else {
          valid_hyperedges[he] = 0;
        }
      }

      // #################### STAGE 3 ####################
      // In the step before we aggregated hyperedges within a bucket data structure.
      // Hyperedges with the same hash/footprint are stored inside the same bucket.
      // We iterate now in parallel over each bucket and sort each bucket
      // after its hash. A bucket is processed by one thread and parallel
      // hyperedges are detected by comparing the pins of hyperedges with
      // the same hash.


      // Helper function that checks if two hyperedges are parallel
      // Note, pins inside the hyperedges are sorted.
      auto check_if_hyperedges_are_parallel = [&](const HyperedgeID lhs,
                                                  const HyperedgeID rhs) {
          const Hyperedge& lhs_he = tmp_hyperedges[lhs];
          const Hyperedge& rhs_he = tmp_hyperedges[rhs];
          if ( lhs_he.size() == rhs_he.size() ) {
            auto lhs_it = lhs_he.pins().begin();
            auto rhs_it = rhs_he.pins().begin();
            for ( ; lhs_it != lhs_he.pins().end(); ++lhs_it, ++rhs_it ) {
              if ( *lhs_it != *rhs_it ) {
                return false;
              }
            }
            return true;
          } else {
            return false;
          }
      };


      for (size_t bucket = 0; bucket < hyperedge_hash_map.numBuckets(); ++bucket) {
        auto &hyperedge_bucket = hyperedge_hash_map.getBucket(bucket);
        std::sort(hyperedge_bucket.begin(), hyperedge_bucket.end(),
                  [&](const ContractedHyperedgeInformation &lhs, const ContractedHyperedgeInformation &rhs) {
                      return std::tie(lhs.hash, lhs.size, lhs.he) < std::tie(rhs.hash, rhs.size, rhs.he);
                  });

        // Parallel Hyperedge Detection
        for (size_t i = 0; i < hyperedge_bucket.size(); ++i) {
          ContractedHyperedgeInformation &contracted_he_lhs = hyperedge_bucket[i];
          if (contracted_he_lhs.valid) {
            const HyperedgeID lhs_he = contracted_he_lhs.he;
            HyperedgeWeight lhs_weight = tmp_hyperedges[lhs_he].weight();
            for (size_t j = i + 1; j < hyperedge_bucket.size(); ++j) {
              ContractedHyperedgeInformation &contracted_he_rhs = hyperedge_bucket[j];
              const HyperedgeID rhs_he = contracted_he_rhs.he;
              if (contracted_he_rhs.valid &&
                  contracted_he_lhs.hash == contracted_he_rhs.hash &&
                  check_if_hyperedges_are_parallel(lhs_he, rhs_he)) {
                // Hyperedges are parallel
                lhs_weight += tmp_hyperedges[rhs_he].weight();
                contracted_he_rhs.valid = false;
                valid_hyperedges[rhs_he] = false;
              } else if (contracted_he_lhs.hash != contracted_he_rhs.hash) {
                // In case, hash of both are not equal we go to the next hyperedge
                // because we compared it with all hyperedges that had an equal hash
                break;
              }
            }
            tmp_hyperedges[lhs_he].setWeight(lhs_weight);
          }
        }
        hyperedge_hash_map.free(bucket);
      }




      MutableHypergraph hypergraph;

      hypergraph._num_hypernodes = unique_communities.size();
      hypergraph._hypernodes.resize(hypergraph._num_hypernodes);
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
      }

      size_t num_valid_hyperedges = 0;
      for (HyperedgeID he : edges()) {
        if (valid_hyperedges[he]) {
          ++num_valid_hyperedges;
        }
      }
      hypergraph._num_hyperedges = num_valid_hyperedges;
      hypergraph._hyperedges.resize(hypergraph._num_hyperedges);
      for (HyperedgeID he : edges()) {
        if (valid_hyperedges[he]) {
          HyperedgeID new_he = hypergraph._hyperedges.size() - num_valid_hyperedges;
          --num_valid_hyperedges;
          hypergraph._hyperedges[new_he] = tmp_hyperedges[he];
          for (const HypernodeID& pin : tmp_hyperedges[he].pins()) {
            hypergraph._hypernodes[pin].addNet(new_he);
            ++hypergraph._num_pins;
            ++hypergraph._total_degree;
          }
          hypergraph._max_edge_size = std::max(static_cast<size_t>(hypergraph._max_edge_size), tmp_hyperedges[he].size());
        }
      }

      // #################### STAGE 4 ####################
      // Set community ids of new nodes

      hypergraph._community_ids.resize(hypergraph._num_hypernodes, 0);
      for (size_t hn : nodes()) {
        if (_community_ids[hn] == static_cast<PartitionID>(kInvalidHypernode) || _hypernodes[hn].is_deleted() || _hypernodes[hn].isDisabled()) {
          continue;
        }
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
        hypergraph._hyperedges[he].enable();
        hypergraph._hyperedges[he].setWeight(_hyperedges[he].weight());
        for (const HypernodeID hn : _hyperedges[he].pins())
        {
          hypergraph._hyperedges[he].addPin(hn);
        }
        if (_hyperedges[he].isDisabled()) {
          hypergraph._hyperedges[he].disable();
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

      hypergraph._community_ids.clear();
      for (size_t hn = 0; hn < _hypernodes.size(); ++hn) {
        if (_hypernodes[hn].is_deleted()) {
          continue;
        }
        hypergraph._community_ids.push_back(_community_ids[hn]);
      }

      //TODO add fixed vertex support

      // create incidence array and incident nets
      hypergraph._incidence_array.resize(_num_pins);
      hypergraph._hyperedges.resize(_num_hyperedges);
      size_t incidence_array_pos = 0;

      // resize mapping vectors
      static_to_mut_he->resize(_hyperedges.size());
      size_t mutable_he_index = -1;
      size_t static_he_index = -1;
      std::vector mut_to_static_he(_hyperedges.size(), kInvalidHyperedge);

      for (size_t he = 0; he < _hyperedges.size(); ++he) {
        ++mutable_he_index;
        if (_hyperedges[he].is_deleted()) {
          deleted_he->push_back(mutable_he_index);
          continue;
        }
        ++static_he_index;
        static_to_mut_he->at(static_he_index) = mutable_he_index;
        mut_to_static_he[mutable_he_index] = static_he_index;
        auto hyperedge_test = StaticHypergraph::Hyperedge();
        hypergraph._hyperedges[static_he_index] = hyperedge_test;
        hypergraph._hyperedges[static_he_index].enable();

        // temporarly enable mutable hypernode to allow copy
        bool is_disabled = !edgeIsEnabled(he);
        if (is_disabled) _hyperedges[he].enable();

        hypergraph._hyperedges[static_he_index].setWeight(_hyperedges[he].weight());
        hypergraph._hyperedges[static_he_index].setSize(_hyperedges[he].size());
        hypergraph._hyperedges[static_he_index].setFirstEntry(incidence_array_pos);

        for (const HypernodeID hn : _hyperedges[he].pins()) {
          hypergraph._incidence_array[incidence_array_pos++] = hn;
        }
        if (is_disabled) {
          hypergraph._hyperedges[static_he_index].disable();
          _hyperedges[he].disable();
        }
      }

      hypergraph._incident_nets.resize(_total_degree);
      hypergraph._hypernodes.resize(_num_hypernodes);
      size_t incident_nets_pos = 0;

      static_to_mut_hn->resize(_hypernodes.size());
      size_t mutable_hn_index = -1;
      size_t static_hn_index = -1;
      std::vector mut_to_static_hn(_hypernodes.size(), kInvalidHypernode);

      for (size_t hn = 0; hn < _hypernodes.size(); ++hn) {
        ++mutable_hn_index;
        if (_hypernodes[hn].is_deleted()) {
          deleted_hn->push_back(mutable_hn_index);
          continue;
        }
        ++static_hn_index;
        static_to_mut_hn->at(static_hn_index) = mutable_hn_index;
        mut_to_static_hn[mutable_hn_index] = static_hn_index;
        hypergraph._hypernodes[static_hn_index] = StaticHypergraph::Hypernode();
        hypergraph._hypernodes[static_hn_index].enable();
        // temporarly enable mutable hypernode to allow copy
        bool is_disabled = !nodeIsEnabled(hn);
        if (is_disabled) _hypernodes[hn].enable();
        hypergraph._hypernodes[static_hn_index].setWeight(_hypernodes[hn].weight());
        hypergraph._hypernodes[static_hn_index].setSize(_hypernodes[hn].size());
        hypergraph._hypernodes[static_hn_index].setFirstEntry(incident_nets_pos);
        for (const HyperedgeID he : _hypernodes[hn].incidentEdges()) {
          hypergraph._incident_nets[incident_nets_pos++] = mut_to_static_he[he];
        }
        if (is_disabled) {
          hypergraph._hypernodes[static_hn_index].disable();
          _hypernodes[hn].disable();
        }
      }

      // map the incidence array to static hypernode ids
      for (size_t i = 0; i < hypergraph._incidence_array.size(); ++i) {
        HypernodeID mut_hn = hypergraph._incidence_array[i];
        hypergraph._incidence_array[i] = mut_to_static_hn[mut_hn];
        ASSERT(hypergraph._incidence_array[i] != kInvalidHypernode);
      }

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

      hypergraph._hypernodes.resize(static_hg._num_hypernodes + deleted_hn->size());

      // TODO reserve instead of single pushes
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

      hypergraph._hyperedges.resize(static_hg._num_hyperedges + deleted_he->size());
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

      return hypergraph;
    }

    MutableHypergraph MutableHypergraph::fromStaticHypergraphSimple(const StaticHypergraph &static_hg) {
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

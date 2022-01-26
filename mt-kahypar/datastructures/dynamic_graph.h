/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Nikolai Maas <nikolai.maas@student.kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <mutex>
#include <queue>

#include "tbb/parallel_for.h"

#include "kahypar/meta/mandatory.h"
#include "kahypar/datastructure/fast_reset_flag_array.h"
#include "kahypar/utils/math.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/dynamic_adjacency_array.h"
#include "mt-kahypar/datastructures/contraction_tree.h"
#include "mt-kahypar/datastructures/thread_safe_fast_reset_flag_array.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/memory_tree.h"

namespace mt_kahypar {
namespace ds {

class DynamicGraph {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  static_assert(std::is_unsigned<HypernodeID>::value, "Hypernode ID must be unsigned");
  static_assert(std::is_unsigned<HyperedgeID>::value, "Hyperedge ID must be unsigned");

  // ! In order to update gain cache correctly for an uncontraction (u,v),
  // ! the partitioned hypergraph has to know wheter v replaces u in a hyperedge
  // ! or both a incident to that hyperedge after uncontraction. Therefore, the partitioned
  // ! hypergraph passes two lambda functions to the batch uncontraction function, one for
  // ! each case.
  using UncontractionFunction = std::function<void (const HypernodeID, const HypernodeID, const HyperedgeID)>;
  #define NOOP_BATCH_FUNC [] (const HypernodeID, const HypernodeID, const HyperedgeID) { }

  // Represents a uncontraction that is assigned to a certain batch
  // and within that batch to a certain position.
  struct BatchAssignment {
    HypernodeID u;
    HypernodeID v;
    size_t batch_index;
    size_t batch_pos;
  };

 private:
  /**
   * Represents a hypernode of the hypergraph and contains all information
   * associated with a vertex.
   */
  class Hypernode {
   public:
    using IDType = HypernodeID;

    Hypernode() :
      _weight(1),
      _community_id(0),
      _batch_idx(std::numeric_limits<HypernodeID>::max()),
      _valid(false) { }

    Hypernode(const bool valid) :
      _weight(1),
      _community_id(0),
      _batch_idx(std::numeric_limits<HypernodeID>::max()),
      _valid(valid) { }

    bool isDisabled() const {
      return _valid == false;
    }

    void enable() {
      ASSERT(isDisabled());
      _valid = true;
    }

    void disable() {
      ASSERT(!isDisabled());
      _valid = false;
    }

    HyperedgeWeight weight() const {
      return _weight;
    }

    void setWeight(HyperedgeWeight weight) {
      ASSERT(!isDisabled());
      _weight = weight;
    }

    PartitionID communityID() const {
      return _community_id;
    }

    void setCommunityID(const PartitionID community_id) {
      ASSERT(!isDisabled());
      _community_id = community_id;
    }

    HypernodeID batchIndex() const {
      return _batch_idx;
    }

    void setBatchIndex(const HypernodeID batch_idx) {
      _batch_idx = batch_idx;
    }

   private:
    // ! Hypernode weight
    HypernodeWeight _weight;
    // ! Community id
    PartitionID _community_id;
    // ! Index of the uncontraction batch in which this hypernode is contained in
    HypernodeID _batch_idx;
    // ! Flag indicating whether or not the element is active.
    bool _valid;
  };

  enum class ContractionResult : uint8_t {
    CONTRACTED = 0,
    PENDING_CONTRACTIONS = 1,
    WEIGHT_LIMIT_REACHED = 2
  };

  using OwnershipVector = parallel::scalable_vector<parallel::IntegralAtomicWrapper<bool>>;
  using ThreadLocalHyperedgeVector = tbb::enumerable_thread_specific<parallel::scalable_vector<HyperedgeID>>;

 public:
  static constexpr bool is_graph = false;
  static constexpr bool is_static_hypergraph = false;
  static constexpr bool is_partitioned = false;

  explicit DynamicGraph() :
    _num_hypernodes(0),
    _num_removed_hypernodes(0),
    _removed_degree_zero_hn_weight(0),
    _num_edges(0),
    _num_removed_hyperedges(0),
    _total_weight(0),
    _version(0),
    _contraction_index(0),
    _hypernodes(),
    _contraction_tree(),
    _acquired_hns(),
    _acquired_hes(),
    _failed_hyperedge_contractions() { }

  DynamicGraph(const DynamicGraph&) = delete;
  DynamicGraph & operator= (const DynamicGraph &) = delete;

  DynamicGraph(DynamicGraph&& other) :
    _num_hypernodes(other._num_hypernodes),
    _num_removed_hypernodes(other._num_removed_hypernodes),
    _removed_degree_zero_hn_weight(other._removed_degree_zero_hn_weight),
    _num_edges(other._num_edges),
    _num_removed_hyperedges(other._num_removed_hyperedges),
    _total_weight(other._total_weight),
    _version(other._version),
    _contraction_index(0),
    _hypernodes(std::move(other._hypernodes)),
    _contraction_tree(std::move(other._contraction_tree)),
    _acquired_hns(std::move(other._acquired_hns)),
    _acquired_hes(std::move(other._acquired_hes)),
    _failed_hyperedge_contractions(std::move(other._failed_hyperedge_contractions)) { }

  DynamicGraph & operator= (DynamicGraph&& other) {
    _num_hypernodes = other._num_hypernodes;
    _num_removed_hypernodes = other._num_removed_hypernodes;
    _num_edges = other._num_edges;
    _num_removed_hyperedges = other._num_removed_hyperedges;
    _removed_degree_zero_hn_weight = other._removed_degree_zero_hn_weight;
    _total_weight = other._total_weight;
    _version = other._version;
    _contraction_index.store(other._contraction_index.load());
    _hypernodes = std::move(other._hypernodes);
    _contraction_tree = std::move(other._contraction_tree);
    _acquired_hns = std::move(other._acquired_hns);
    _acquired_hes = std::move(other._acquired_hes);
    _failed_hyperedge_contractions = std::move(other._failed_hyperedge_contractions);
    return *this;
  }

  ~DynamicGraph() {
    freeInternalData();
  }

  // ####################### General Hypergraph Stats #######################

  // ! Initial number of hypernodes
  HypernodeID initialNumNodes() const {
    return _num_hypernodes;
  }

  // ! Number of removed hypernodes
  HypernodeID numRemovedHypernodes() const {
    return _num_removed_hypernodes;
  }

  // ! Weight of removed degree zero vertics
  HypernodeWeight weightOfRemovedDegreeZeroVertices() const {
    return _removed_degree_zero_hn_weight;
  }

  // ! Initial number of hyperedges
  HyperedgeID initialNumEdges() const {
    return _num_edges;
  }

  // ! Number of removed hyperedges
  HyperedgeID numRemovedHyperedges() const {
    return _num_removed_hyperedges;
  }

  // ! Set the number of removed hyperedges
  void setNumRemovedHyperedges(const HyperedgeID num_removed_hyperedges) {
    _num_removed_hyperedges = num_removed_hyperedges;
  }

  // ! Initial number of pins
  HypernodeID initialNumPins() const {
    return _num_edges;
  }

  // ! Initial sum of the degree of all vertices
  HypernodeID initialTotalVertexDegree() const {
    return _num_edges;
  }

  // ! Total weight of hypergraph
  HypernodeWeight totalWeight() const {
    return _total_weight;
  }

  // ! Recomputes the total weight of the hypergraph (parallel)
  void updateTotalWeight(parallel_tag_t);

  // ! Recomputes the total weight of the hypergraph (sequential)
  void updateTotalWeight();

  // ####################### Iterators #######################

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const F& f) {
    static_cast<const DynamicGraph&>(*this).doParallelForAllNodes(f);
  }

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const F& f) const {
    tbb::parallel_for(ID(0), _num_hypernodes, [&](const HypernodeID& hn) {
      if ( nodeIsEnabled(hn) ) {
        f(hn);
      }
    });
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const F& f) {
    static_cast<const DynamicGraph&>(*this).doParallelForAllEdges(f);
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const F& f) const {
    tbb::parallel_for(ID(0), _num_edges, [&](const HyperedgeID& he) {
      if ( edgeIsEnabled(he) ) {
        f(he);
      }
    });
  }

  // ####################### Hypernode Information #######################

  // ! Weight of a vertex
  HypernodeWeight nodeWeight(const HypernodeID u) const {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return hypernode(u).weight();
  }

  // ! Sets the weight of a vertex
  void setNodeWeight(const HypernodeID u, const HypernodeWeight weight) {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).setWeight(weight);
  }

  // ! Degree of a hypernode
  HyperedgeID nodeDegree(const HypernodeID u) const {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return _incident_edges.nodeDegree(u);
  }

  // ! Returns, whether a hypernode is enabled or not
  bool nodeIsEnabled(const HypernodeID u) const {
    return !hypernode(u).isDisabled();
  }

  // ! Enables a hypernode (must be disabled before)
  void enableHypernode(const HypernodeID u) {
    hypernode(u).enable();
  }

  // ! Disables a hypernode (must be enabled before)
  void disableHypernode(const HypernodeID u) {
    hypernode(u).disable();
  }

  // ! Removes a hypernode (must be enabled before)
  void removeHypernode(const HypernodeID u) {
    hypernode(u).disable();
    ++_num_removed_hypernodes;
  }

  // ! Removes a degree zero hypernode
  void removeDegreeZeroHypernode(const HypernodeID u) {
    ASSERT(nodeDegree(u) == 0);
    removeHypernode(u);
    _removed_degree_zero_hn_weight += nodeWeight(u);
  }

  // ! Restores a degree zero hypernode
  void restoreDegreeZeroHypernode(const HypernodeID u) {
    hypernode(u).enable();
    ASSERT(nodeDegree(u) == 0);
    _removed_degree_zero_hn_weight -= nodeWeight(u);
  }

  // ####################### Hyperedge Information #######################

  // ! Weight of a hyperedge
  // HypernodeWeight edgeWeight(const HyperedgeID e) const {
  //   ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
  //   return hyperedge(e).weight();
  // }

  // ! Sets the weight of a hyperedge
  // void setEdgeWeight(const HyperedgeID e, const HyperedgeWeight weight) {
  //   ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
  //   return hyperedge(e).setWeight(weight);
  // }

  // ! Number of pins of a hyperedge
  HypernodeID edgeSize(const HyperedgeID e) const {
    return 2;
  }

  // ! Maximum size of a hyperedge
  HypernodeID maxEdgeSize() const {
    return 2;
  }

  // ! Hash value defined over the pins of a hyperedge
  // size_t edgeHash(const HyperedgeID e) const {
  //   ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
  //   return hyperedge(e).hash();
  // }

  // ! Returns, whether a hyperedge is enabled or not
  // bool edgeIsEnabled(const HyperedgeID e) const {
  //   return !hyperedge(e).isDisabled();
  // }

  // ! Enables a hyperedge (must be disabled before)
  // void enableHyperedge(const HyperedgeID e) {
  //   hyperedge(e).enable();
  // }

  // ! Disabled a hyperedge (must be enabled before)
  // void disableHyperedge(const HyperedgeID e) {
  //   hyperedge(e).disable();
  // }

  // ####################### Community Information #######################

  // ! Community id which hypernode u is assigned to
  PartitionID communityID(const HypernodeID u) const {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return hypernode(u).communityID();
  }

  // ! Assign a community to a hypernode
  // ! Note, in order to use all community-related functions, initializeCommunities()
  // ! have to be called after assigning to each vertex a community id
  void setCommunityID(const HypernodeID u, const PartitionID community_id) {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).setCommunityID(community_id);
  }

  // ####################### Contract / Uncontract #######################

  DynamicGraph contract(parallel::scalable_vector<HypernodeID>&) {
    ERROR("contract(c, id) is not supported in dynamic graph");
    return DynamicGraph();
  }

  /**!
   * Registers a contraction in the hypergraph whereas vertex u is the representative
   * of the contraction and v its contraction partner. Several threads can call this function
   * in parallel. The function adds the contraction of u and v to a contraction tree that determines
   * a parallel execution order and synchronization points for all running contractions.
   * The contraction can be executed by calling function contract(v, max_node_weight).
   */
  bool registerContraction(const HypernodeID u, const HypernodeID v);

  /**!
   * Contracts a previously registered contraction. Representative u of vertex v is looked up
   * in the contraction tree and performed if there are no pending contractions in the subtree
   * of v and the contractions respects the maximum allowed node weight. If (u,v) is the last
   * pending contraction in the subtree of u then the function recursively contracts also
   * u (if any contraction is registered). Therefore, function can return several contractions
   * or also return an empty contraction vector.
   */
  size_t contract(const HypernodeID v,
                  const HypernodeWeight max_node_weight = std::numeric_limits<HypernodeWeight>::max());

  /**
   * Uncontracts a batch of contractions in parallel. The batches must be uncontracted exactly
   * in the order computed by the function createBatchUncontractionHierarchy(...).
   * The two uncontraction functions are required by the partitioned hypergraph to restore
   * pin counts and gain cache values.
   */
  void uncontract(const Batch& batch,
                  const UncontractionFunction& case_one_func = NOOP_BATCH_FUNC,
                  const UncontractionFunction& case_two_func = NOOP_BATCH_FUNC);

  /**
   * Computes a batch uncontraction hierarchy. A batch is a vector of mementos
   * (uncontractions) that are uncontracted in parallel. The function returns a vector
   * of versioned batch vectors. A new version of the hypergraph is induced if we perform
   * single-pin and parallel net detection. Once we process all batches of a versioned
   * batch vector, we have to restore all previously removed single-pin and parallel nets
   * in order to uncontract the next batch vector. We create for each version of the
   * hypergraph a seperate batch uncontraction hierarchy (see createBatchUncontractionHierarchyOfVersion(...))
   */
  VersionedBatchVector createBatchUncontractionHierarchy(const size_t batch_size,
                                                         const bool test = false);

  // ! Only for testing
  VersionedBatchVector createBatchUncontractionHierarchy(ContractionTree&& tree,
                                                         const size_t batch_size,
                                                         const size_t num_versions = 1) {
    ASSERT(num_versions > 0);
    _version = num_versions - 1;
    _contraction_tree = std::move(tree);
    return createBatchUncontractionHierarchy(batch_size, true);
  }

  // ! Only for testing
  HypernodeID contractionTree(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return _contraction_tree.parent(u);
  }

  // ! Only for testing
  HypernodeID pendingContractions(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return _contraction_tree.pendingContractions(u);
  }

  // ! Only for testing
  void decrementPendingContractions(const HypernodeID u) {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    _contraction_tree.decrementPendingContractions(u);
  }

  // ####################### Remove / Restore Hyperedges #######################

  /*!
  * Removes a hyperedge from the hypergraph. This includes the removal of he from all
  * of its pins and to disable the hyperedge.
  *
  * NOTE, this function is not thread-safe and should only be called in a single-threaded
  * setting.
  */
  // void removeEdge(const HyperedgeID he) {
  //   ASSERT(edgeIsEnabled(he), "Hyperedge" << he << "is disabled");
  //   kahypar::ds::FastResetFlagArray<>& he_to_remove = _he_bitset.local();
  //   he_to_remove.set(he, true);
  //   for ( const HypernodeID& pin : pins(he) ) {
  //     _incident_nets.removeIncidentNets(pin, he_to_remove);
  //   }
  //   ++_num_removed_hyperedges;
  //   disableHyperedge(he);
  // }

  /*!
  * Removes a hyperedge from the hypergraph. This includes the removal of he from all
  * of its pins and to disable the hyperedge. Note, in contrast to removeEdge, this function
  * removes hyperedge from all its pins in parallel.
  *
  * NOTE, this function is not thread-safe and should only be called in a single-threaded
  * setting.
  */
  void removeLargeEdge(const HyperedgeID he) {
    ERROR("removeLargeEdge is not supported in dynamic graph");
  }

  /*!
   * Restores a large hyperedge previously removed from the hypergraph.
   */
  void restoreLargeEdge(const HyperedgeID& he) {
    ERROR("restoreLargeEdge is not supported in dynamic graph");
  }

  /**
   * Removes single-pin and parallel nets from the hypergraph. The weight
   * of a set of identical nets is aggregated in one representative hyperedge
   * and single-pin hyperedges are removed. Returns a vector of removed hyperedges.
   */
  parallel::scalable_vector<ParallelHyperedge> removeSinglePinAndParallelHyperedges();

  /**
   * Restores a previously removed set of singple-pin and parallel hyperedges. Note, that hes_to_restore
   * must be exactly the same and given in the reverse order as returned by removeSinglePinAndParallelNets(...).
   */
  void restoreSinglePinAndParallelNets(const parallel::scalable_vector<ParallelHyperedge>& hes_to_restore);

  // ####################### Initialization / Reset Functions #######################

  // ! Reset internal community information
  void setCommunityIDs(const parallel::scalable_vector<PartitionID>& community_ids) {
    ASSERT(community_ids.size() == UI64(_num_hypernodes));
    doParallelForAllNodes([&](const HypernodeID& hn) {
      hypernode(hn).setCommunityID(community_ids[hn]);
    });

  }

  // ####################### Copy #######################

  // ! Copy dynamic hypergraph in parallel
  DynamicGraph copy(parallel_tag_t);

  // ! Copy dynamic hypergraph sequential
  DynamicGraph copy();

  // ! Reset internal data structure
  void reset() {
    _contraction_tree.reset();
    _incident_edges.reset();
    _version = 0;
  }

  // ! Free internal data in parallel
  void freeInternalData() {
    _num_hypernodes = 0;
    _num_edges = 0;
  }

  void freeTmpContractionBuffer() {
    ERROR("freeTmpContractionBuffer() is not supported in dynamic hypergraph");
  }

  void memoryConsumption(utils::MemoryTreeNode* parent) const;

  // ! Only for testing
  bool verifyIncidenceArrayAndIncidentNets();

 private:
  // friend class DynamicGraphFactory;
  template<typename Hypergraph>
  friend class CommunitySupport;
  template <typename Hypergraph,
            typename HypergraphFactory>
  friend class PartitionedHypergraph;

  // ####################### Acquiring / Releasing Ownership #######################

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void acquireHypernode(const HypernodeID u) {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    bool expected = false;
    bool desired = true;
    while ( !_acquired_hns[u].compare_exchange_strong(expected, desired) ) {
      expected = false;
    }
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool tryAcquireHypernode(const HypernodeID u) {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    bool expected = false;
    bool desired = true;
    return _acquired_hns[u].compare_exchange_strong(expected, desired);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void releaseHypernode(const HypernodeID u) {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    ASSERT(_acquired_hns[u], "Hypernode" << u << "is not acquired!");
    _acquired_hns[u] = false;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void acquireHyperedge(const HyperedgeID e) {
    ASSERT(e < _num_edges, "Hyperedge" << e << "does not exist");
    bool expected = false;
    bool desired = true;
    while ( !_acquired_hes[e].compare_exchange_strong(expected, desired) ) {
      expected = false;
    }
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool tryAcquireHyperedge(const HyperedgeID e) {
    ASSERT(e < _num_edges, "Hyperedge" << e << "does not exist");
    bool expected = false;
    bool desired = true;
    return _acquired_hes[e].compare_exchange_strong(expected, desired);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void releaseHyperedge(const HyperedgeID e) {
    ASSERT(e < _num_edges, "Hyperedge" << e << "does not exist");
    ASSERT(_acquired_hes[e], "Hyperedge" << e << "is not acquired!");
    _acquired_hes[e] = false;
  }

  // ####################### Hypernode Information #######################

  // ! Accessor for hypernode-related information
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Hypernode& hypernode(const HypernodeID u) const {
    ASSERT(u <= _num_hypernodes, "Hypernode" << u << "does not exist");
    return _hypernodes[u];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hypernode& hypernode(const HypernodeID u) {
    return const_cast<Hypernode&>(static_cast<const DynamicGraph&>(*this).hypernode(u));
  }

  // MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE IteratorRange<IncidentNetsIterator> incident_nets_of(const HypernodeID u,
  //                                                                                         const size_t pos = 0) const {
  //   return _incident_nets.incidentEdges(u, pos);
  // }

  // ####################### Hyperedge Information #######################

  // ! Accessor for hyperedge-related information
  // MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Hyperedge& hyperedge(const HyperedgeID e) const {
  //   ASSERT(e <= _num_edges, "Hyperedge" << e << "does not exist");
  //   return _hyperedges[e];
  // }

  // ! To avoid code duplication we implement non-const version in terms of const version
  // MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hyperedge& hyperedge(const HyperedgeID e) {
  //   return const_cast<Hyperedge&>(static_cast<const DynamicGraph&>(*this).hyperedge(e));
  // }

  // ####################### Contract / Uncontract #######################

  /**!
   * Contracts a previously registered contraction. The contraction of u and v is
   * performed if there are no pending contractions in the subtree of v and the
   * contractions respects the maximum allowed node weight. In case the contraction
   * was performed successfully, enum type CONTRACTED is returned. If contraction
   * was not performed either WEIGHT_LIMIT_REACHED (in case sum of both vertices is
   * greater than the maximum allowed node weight) or PENDING_CONTRACTIONS (in case
   * there are some unfinished contractions in the subtree of v) is returned.
   */
  ContractionResult contract(const HypernodeID u,
                             const HypernodeID v,
                             const HypernodeWeight max_node_weight);

  // ! Performs the contraction of (u,v) inside hyperedge he
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void contractHyperedge(const HypernodeID u, const HypernodeID v, const HyperedgeID he,
                                                            kahypar::ds::FastResetFlagArray<>& shared_incident_nets_u_and_v);

  // ! Restore the size of the hyperedge to the size before the batch with
  // ! index batch_index was contracted.
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void restoreHyperedgeSizeForBatch(const HyperedgeID he,
                                                                       const HypernodeID batch_index);

  // ! Search for the position of pin u in hyperedge he in the incidence array
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE size_t findPositionOfPinInIncidenceArray(const HypernodeID u,
                                                                              const HyperedgeID he);

  // bool verifyBatchIndexAssignments(
  //   const BatchIndexAssigner& batch_assigner,
  //   const parallel::scalable_vector<parallel::scalable_vector<BatchAssignment>>& local_batch_assignments) const;

  /**
   * Computes a batch uncontraction hierarchy for a specific version of the hypergraph.
   * A batch is a vector of mementos (uncontractions) that are uncontracted in parallel.
   * Each time we perform single-pin and parallel net detection we create a new version of
   * the hypergraph.
   * A batch of uncontractions that is uncontracted in parallel must satisfy two conditions:
   *  1.) All representatives must be active vertices of the hypergraph
   *  2.) For a specific representative its contraction partners must be uncontracted in reverse
   *      contraction order. Meaning that a contraction (u, v) that happens before a contraction (u, w)
   *      must be uncontracted in a batch that is part of the same batch or a batch uncontracted after the
   *      batch which (u, w) is part of. This ensures that a parallel batch uncontraction does not
   *      increase the objective function.
   * We use the contraction tree to create a batch uncontraction order. Note, uncontractions from
   * different subtrees can be interleaved abitrary. To ensure condition 1.) we peform a BFS starting
   * from all roots of the contraction tree. Each BFS level induces a new batch. Since we contract
   * vertices in parallel its not possible to create a relative order of the contractions which is
   * neccassary for condition 2.). However, during a contraction we store a start and end "timepoint"
   * of a contraction. If two contractions time intervals do not intersect, we can determine
   * which contraction happens strictly before the other. If they intersect, it is not possible to
   * give a relative order. To ensure condition 2.) we sort the childs of a vertex in the contraction tree
   * after its time intervals. Once we add a uncontraction (u,v) to a batch, we also add all uncontractions
   * (u,w) to the batch which intersect the time interval of (u,v). To merge uncontractions of different
   * subtrees in a batch, we insert all eligble uncontractions into a max priority queue with the subtree
   * size of the contraction partner as key. We insert uncontractions into the current batch as long
   * as the maximum batch size is not reached or the PQ is empty. Once the batch reaches its maximum
   * batch size, we create a new empty batch. If the PQ is empty, we replace it with the PQ of the next
   * BFS level. With this approach heavy vertices are uncontracted earlier (subtree size in the PQ as key = weight of
   * a vertex for an unweighted hypergraph) such that average node weight of the hypergraph decreases faster and
   * local searches are more effective in early stages of the uncontraction hierarchy where hyperedge sizes are
   * usually smaller than on the original hypergraph.
   */
  // BatchVector createBatchUncontractionHierarchyForVersion(BatchIndexAssigner& batch_assigner,
  //                                                         const size_t version);

  // ! Number of hypernodes
  HypernodeID _num_hypernodes;
  // ! Number of removed hypernodes
  HypernodeID _num_removed_hypernodes;
  // ! Number of removed degree zero hypernodes
  HypernodeWeight _removed_degree_zero_hn_weight;
  // ! Number of hyperedges
  HyperedgeID _num_edges;
  // ! Number of removed hyperedges
  HyperedgeID _num_removed_hyperedges;
  // ! Total weight of hypergraph
  HypernodeWeight _total_weight;
  // ! Version of the hypergraph, each time we remove a single-pin and parallel nets,
  // ! we create a new version
  size_t _version;
  // ! Contraction Index, increment whenever a contraction terminates
  std::atomic<HypernodeID> _contraction_index;

  // ! Hypernodes
  Array<Hypernode> _hypernodes;
  // ! Contraction Tree
  ContractionTree _contraction_tree;
  // ! Pins of hyperedges
  IncidentEdgeArray _incident_edges;
  // ! Atomic bool vector used to acquire unique ownership of hypernodes
  OwnershipVector _acquired_hns;

  // ! Atomic bool vector used to acquire unique ownership of hyperedges
  OwnershipVector _acquired_hes;
  // ! Collects hyperedge contractions that failed due to failed acquired ownership
  ThreadLocalHyperedgeVector _failed_hyperedge_contractions;
};

} // namespace ds
} // namespace mt_kahypar
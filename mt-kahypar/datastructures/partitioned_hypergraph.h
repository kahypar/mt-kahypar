/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <atomic>
#include <type_traits>
#include <mutex>

#include "tbb/parallel_invoke.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/connectivity_set.h"
#include "mt-kahypar/datastructures/pin_count_in_part.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/stl/thread_locals.h"
#include "mt-kahypar/utils/range.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
namespace ds {

template <typename Hypergraph = Mandatory,
          typename HypergraphFactory = Mandatory>
class PartitionedHypergraph {
private:
  static_assert(!Hypergraph::is_partitioned,  "Only unpartitioned hypergraphs are allowed");

  // ! Function that will be called for each incident hyperedge of a moved vertex with the following arguments
  // !  1) hyperedge ID, 2) weight, 3) size, 4) pin count in from-block after move, 5) pin count in to-block after move
  // ! Can be implemented to obtain correct km1 or cut improvements of the move
  using DeltaFunction = std::function<void (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID)>;
  #define NOOP_FUNC [] (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID) { }

  // REVIEW NOTE: Can't we use a lambda in changeNodePart. And write a second function that calls the first with a lambda that does nothing.
  // Then we could guarantee inlining
  // This would also reduce the code/documentation copy-pasta for with or without gain updates

  static constexpr bool enable_heavy_assert = false;

 public:
  static constexpr bool is_static_hypergraph = Hypergraph::is_static_hypergraph;
  static constexpr bool is_partitioned = true;
  static constexpr bool supports_connectivity_set = true;

  static constexpr HyperedgeID HIGH_DEGREE_THRESHOLD = ID(100000);

  using HypernodeIterator = typename Hypergraph::HypernodeIterator;
  using HyperedgeIterator = typename Hypergraph::HyperedgeIterator;
  using IncidenceIterator = typename Hypergraph::IncidenceIterator;
  using IncidentNetsIterator = typename Hypergraph::IncidentNetsIterator;

  PartitionedHypergraph() = default;

  explicit PartitionedHypergraph(const PartitionID k,
                                 Hypergraph& hypergraph) :
    _is_gain_cache_initialized(false),
    _top_level_num_nodes(hypergraph.initialNumNodes()),
    _k(k),
    _hg(&hypergraph),
    _part_weights(k, CAtomic<HypernodeWeight>(0)),
    _part_ids(
        "Refinement", "part_ids", hypergraph.initialNumNodes(), false, false),
    _pins_in_part(hypergraph.initialNumEdges(), k, hypergraph.maxEdgeSize(), false),
    _connectivity_set(hypergraph.initialNumEdges(), k, false),
    _gain_cache(),
    _pin_count_update_ownership(
        "Refinement", "pin_count_update_ownership", hypergraph.initialNumEdges(), true, false) {
    _part_ids.assign(hypergraph.initialNumNodes(), kInvalidPartition, false);
  }

  explicit PartitionedHypergraph(const PartitionID k,
                                 Hypergraph& hypergraph,
                                 parallel_tag_t) :
    _is_gain_cache_initialized(false),
    _top_level_num_nodes(hypergraph.initialNumNodes()),
    _k(k),
    _hg(&hypergraph),
    _part_weights(k, CAtomic<HypernodeWeight>(0)),
    _part_ids(),
    _pins_in_part(),
    _connectivity_set(0, 0),
    _gain_cache(),
    _pin_count_update_ownership() {
    tbb::parallel_invoke([&] {
      _part_ids.resize(
        "Refinement", "vertex_part_info", hypergraph.initialNumNodes());
      _part_ids.assign(hypergraph.initialNumNodes(), kInvalidPartition);
    }, [&] {
      _pins_in_part.initialize(hypergraph.initialNumEdges(), k, hypergraph.maxEdgeSize());
    }, [&] {
      _connectivity_set = ConnectivitySets(hypergraph.initialNumEdges(), k);
    }, [&] {
      _pin_count_update_ownership.resize(
        "Refinement", "pin_count_update_ownership", hypergraph.initialNumEdges(), true);
    });
  }

  // REVIEW NOTE why do we delete copy assignment/construction? wouldn't it be useful to make a copy, e.g. for initial partitioning
  PartitionedHypergraph(const PartitionedHypergraph&) = delete;
  PartitionedHypergraph & operator= (const PartitionedHypergraph &) = delete;

  PartitionedHypergraph(PartitionedHypergraph&& other) = default;
  PartitionedHypergraph & operator= (PartitionedHypergraph&& other) = default;

  ~PartitionedHypergraph() {
    freeInternalData();
  }

  void resetData() {
    _is_gain_cache_initialized = false;
    tbb::parallel_invoke([&] {
    }, [&] {
      _part_ids.assign(_part_ids.size(), kInvalidPartition);
    }, [&] {
      _pins_in_part.data().assign(_pins_in_part.data().size(), 0);
    }, [&] {
      _connectivity_set.reset();
    }, [&] {
      for (auto& x : _part_weights) x.store(0, std::memory_order_relaxed);
    });
  }

  // ####################### General Hypergraph Stats ######################

  Hypergraph& hypergraph() {
    ASSERT(_hg);
    return *_hg;
  }

  void setHypergraph(Hypergraph& hypergraph) {
    _hg = &hypergraph;
  }

  // ! Initial number of hypernodes
  HypernodeID initialNumNodes() const {
    return _hg->initialNumNodes();
  }

  // ! Number of removed hypernodes
  HypernodeID numRemovedHypernodes() const {
    return _hg->numRemovedHypernodes();
  }

  // ! Initial number of hyperedges
  HyperedgeID initialNumEdges() const {
    return _hg->initialNumEdges();
  }

  // ! Initial number of pins
  HypernodeID initialNumPins() const {
    return _hg->initialNumPins();
  }

  // ! Initial sum of the degree of all vertices
  HypernodeID initialTotalVertexDegree() const {
    return _hg->initialTotalVertexDegree();
  }

  // ! Total weight of hypergraph
  HypernodeWeight totalWeight() const {
    return _hg->totalWeight();
  }

  // ! Number of blocks this hypergraph is partitioned into
  PartitionID k() const {
    return _k;
  }


  // ####################### Iterators #######################

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const F& f) {
    static_cast<const PartitionedHypergraph&>(*this).doParallelForAllNodes(f);
  }

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const F& f) const {
    _hg->doParallelForAllNodes(f);
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const F& f) {
    static_cast<const PartitionedHypergraph&>(*this).doParallelForAllEdges(f);
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const F& f) const {
    _hg->doParallelForAllEdges(f);
  }

  // ! Returns an iterator over the set of active nodes of the hypergraph
  IteratorRange<HypernodeIterator> nodes() const {
    return _hg->nodes();
  }

  // ! Returns an iterator over the set of active edges of the hypergraph
  IteratorRange<HyperedgeIterator> edges() const {
    return _hg->edges();
  }

  // ! Returns a range to loop over the incident nets of hypernode u.
  IteratorRange<IncidentNetsIterator> incidentEdges(const HypernodeID u) const {
    return _hg->incidentEdges(u);
  }

  // ! Returns a range to loop over the pins of hyperedge e.
  IteratorRange<IncidenceIterator> pins(const HyperedgeID e) const {
    return _hg->pins(e);
  }

  // ! Returns a range to loop over the set of block ids contained in hyperedge e.
  IteratorRange<ConnectivitySets::Iterator> connectivitySet(const HyperedgeID e) const {
    ASSERT(_hg->edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    return _connectivity_set.connectivitySet(e);
  }

  // ####################### Hypernode Information #######################

  // ! Weight of a vertex
  HypernodeWeight nodeWeight(const HypernodeID u) const {
    return _hg->nodeWeight(u);
  }

  // ! Sets the weight of a vertex
  void setNodeWeight(const HypernodeID u, const HypernodeWeight weight) {
    const PartitionID block = partID(u);
    if ( block != kInvalidPartition ) {
      ASSERT(block < _k);
      const HypernodeWeight delta = weight - _hg->nodeWeight(u);
      _part_weights[block] += delta;
    }
    _hg->setNodeWeight(u, weight);
  }

  // ! Degree of a hypernode
  HyperedgeID nodeDegree(const HypernodeID u) const {
    return _hg->nodeDegree(u);
  }

  // ! Returns, whether a hypernode is enabled or not
  bool nodeIsEnabled(const HypernodeID u) const {
    return _hg->nodeIsEnabled(u);
  }

  // ! Enables a hypernode (must be disabled before)
  void enableHypernode(const HypernodeID u) {
    _hg->enableHypernode(u);
  }

  // ! Disable a hypernode (must be enabled before)
  void disableHypernode(const HypernodeID u) {
    _hg->disableHypernode(u);
  }

  // ! Restores a degree zero hypernode
  void restoreDegreeZeroHypernode(const HypernodeID u, const PartitionID to) {
    _hg->restoreDegreeZeroHypernode(u);
    setNodePart(u, to);
  }

  // ####################### Hyperedge Information #######################

  // ! Weight of a hyperedge
  HypernodeWeight edgeWeight(const HyperedgeID e) const {
    return _hg->edgeWeight(e);
  }

  // ! Unique id of a hyperedge
  HyperedgeID uniqueEdgeID(const HyperedgeID e) const {
    return e;
  }

  // ! Sets the weight of a hyperedge
  void setEdgeWeight(const HyperedgeID e, const HyperedgeWeight weight) {
    _hg->setEdgeWeight(e, weight);
  }

  // ! Number of pins of a hyperedge
  HypernodeID edgeSize(const HyperedgeID e) const {
    return _hg->edgeSize(e);
  }

  // ! Returns, whether a hyperedge is enabled or not
  bool edgeIsEnabled(const HyperedgeID e) const {
    return _hg->edgeIsEnabled(e);
  }

  // ! Enables a hyperedge (must be disabled before)
  void enableHyperedge(const HyperedgeID e) {
    _hg->enableHyperedge(e);
  }

  // ! Disabled a hyperedge (must be enabled before)
  void disableHyperedge(const HyperedgeID e) {
    _hg->disableHyperedge(e);
  }

  // ####################### Uncontraction #######################

  /**
   * Uncontracts a batch of contractions in parallel. The batches must be uncontracted exactly
   * in the order computed by the function createBatchUncontractionHierarchy(...).
   */
  void uncontract(const Batch& batch) {
    // Set block ids of contraction partners
    tbb::parallel_for(0UL, batch.size(), [&](const size_t i) {
      const Memento& memento = batch[i];
      ASSERT(nodeIsEnabled(memento.u));
      ASSERT(!nodeIsEnabled(memento.v));
      const PartitionID part_id = partID(memento.u);
      ASSERT(part_id != kInvalidPartition && part_id < _k);
      setOnlyNodePart(memento.v, part_id);
    });

    _hg->uncontract(batch,
      [&](const HypernodeID u, const HypernodeID v, const HyperedgeID he) {
        // In this case, u and v are incident to hyperedge he after uncontraction
        const PartitionID block = partID(u);
        const HypernodeID pin_count_in_part_after = incrementPinCountInPartWithoutGainUpdate(he, block);
        ASSERT(pin_count_in_part_after > 1, V(u) << V(v) << V(he));

        if ( _is_gain_cache_initialized ) {
          // If u was the only pin of hyperedge he in its block before then moving out vertex u
          // of hyperedge he does not decrease the connectivity any more after the
          // uncontraction => b(u) -= w(he)
          const HyperedgeWeight edge_weight = edgeWeight(he);
          if ( pin_count_in_part_after == 2 ) {
            // u might be replaced by an other vertex in the batch
            // => search for other pin of the corresponding block and
            // substract edge weight.
            for ( const HypernodeID& pin : pins(he) ) {
              if ( pin != v && partID(pin) == block ) {
                _gain_cache[penalty_index(pin)].add_fetch(edge_weight, std::memory_order_relaxed);
                break;
              }
            }
          }

          _gain_cache[penalty_index(v)].add_fetch(edge_weight, std::memory_order_relaxed);
          // For all blocks contained in the connectivity set of hyperedge he
          // we increase the move_to_benefit for vertex v by w(e)
          for ( const PartitionID block : _connectivity_set.connectivitySet(he) ) {
            _gain_cache[benefit_index(v, block)].add_fetch(
              edge_weight, std::memory_order_relaxed);
          }
        }
      },
      [&](const HypernodeID u, const HypernodeID v, const HyperedgeID he) {
        // In this case, u is replaced by v in hyperedge he
        // => Pin counts of hyperedge he does not change
        if ( _is_gain_cache_initialized ) {
          const PartitionID block = partID(u);
          const HyperedgeWeight edge_weight = edgeWeight(he);
          // Since u is no longer incident to hyperedge he its contribution for decreasing
          // the connectivity of he is shifted to vertex v
          if ( pinCountInPart(he, block) == 1 ) {
            _gain_cache[penalty_index(u)].add_fetch(edge_weight, std::memory_order_relaxed);
            _gain_cache[penalty_index(v)].sub_fetch(edge_weight, std::memory_order_relaxed);
          }

          _gain_cache[penalty_index(u)].sub_fetch(
            edge_weight, std::memory_order_relaxed);
          _gain_cache[penalty_index(v)].add_fetch(
            edge_weight, std::memory_order_relaxed);
          // For all blocks contained in the connectivity set of hyperedge he
          // we increase the move_to_benefit for vertex v by w(e) and decrease
          // it for vertex u by w(e)
          for ( const PartitionID block : _connectivity_set.connectivitySet(he) ) {
            _gain_cache[benefit_index(u, block)].sub_fetch(
              edge_weight, std::memory_order_relaxed);
            _gain_cache[benefit_index(v, block)].add_fetch(
              edge_weight, std::memory_order_relaxed);
          }
        }
      });
  }

  // ####################### Restore Hyperedges #######################

  /*!
   * Restores a large hyperedge previously removed from the hypergraph.
   */
  void restoreLargeEdge(const HyperedgeID& he) {
    _hg->restoreLargeEdge(he);

    // Recalculate pin count in parts
    const size_t incidence_array_start = _hg->hyperedge(he).firstEntry();
    const size_t incidence_array_end = _hg->hyperedge(he).firstInvalidEntry();
    tls_enumerable_thread_specific< vec<HypernodeID> > ets_pin_count_in_part(_k, 0);
    tbb::parallel_for(incidence_array_start, incidence_array_end, [&](const size_t pos) {
      const HypernodeID pin = _hg->_incidence_array[pos];
      const PartitionID block = partID(pin);
      ++ets_pin_count_in_part.local()[block];
    });

    // Aggregate local pin count for each block
    for ( PartitionID block = 0; block < _k; ++block ) {
      HypernodeID pin_count_in_part = 0;
      for ( const vec<HypernodeID>& local_pin_count : ets_pin_count_in_part ) {
        pin_count_in_part += local_pin_count[block];
      }

      if ( pin_count_in_part > 0 ) {
        _pins_in_part.setPinCountInPart(he, block, pin_count_in_part);
        _connectivity_set.add(he, block);
      }
    }
  }

  /**
   * Restores a previously removed set of singple-pin and parallel hyperedges. Note, that hes_to_restore
   * must be exactly the same and given in the reverse order as returned by removeSinglePinAndParallelNets(...).
   */
  void restoreSinglePinAndParallelNets(const parallel::scalable_vector<typename Hypergraph::ParallelHyperedge>& hes_to_restore) {
    // Restore hyperedges in hypergraph
    _hg->restoreSinglePinAndParallelNets(hes_to_restore);

    // Compute pin counts of restored hyperedges and gain cache values of vertices contained
    // single-pin hyperedges. Note, that restoring parallel hyperedges does not change any
    // value in the gain cache, since it already contributes to the gain via its representative.
    tls_enumerable_thread_specific< vec<HypernodeID> > ets_pin_count_in_part(_k, 0);
    tbb::parallel_for(0UL, hes_to_restore.size(), [&](const size_t i) {
      const HyperedgeID he = hes_to_restore[i].removed_hyperedge;
      const HyperedgeID representative = hes_to_restore[i].representative;
      ASSERT(edgeIsEnabled(he));
      const bool is_single_pin_he = edgeSize(he) == 1;
      if ( is_single_pin_he ) {
        // Restore single-pin net
        HypernodeID single_vertex_of_he = kInvalidHypernode;
        for ( const HypernodeID& pin : pins(he) ) {
          single_vertex_of_he = pin;
        }
        ASSERT(single_vertex_of_he != kInvalidHypernode);

        const PartitionID block_of_single_pin = partID(single_vertex_of_he);
        _connectivity_set.add(he, block_of_single_pin);
        _pins_in_part.setPinCountInPart(he, block_of_single_pin, 1);

        if ( _is_gain_cache_initialized ) {
          _gain_cache[benefit_index(
            single_vertex_of_he, block_of_single_pin)].add_fetch(
              edgeWeight(he), std::memory_order_relaxed);
        }
      } else {
        // Restore parallel net => pin count information given by representative
        ASSERT(edgeIsEnabled(representative));
        for ( const PartitionID& block : connectivitySet(representative) ) {
          _connectivity_set.add(he, block);
          _pins_in_part.setPinCountInPart(he, block, pinCountInPart(representative, block));
        }

        HEAVY_REFINEMENT_ASSERT([&] {
          for ( PartitionID block = 0; block < _k; ++block ) {
            if ( pinCountInPart(he, block) != pinCountInPartRecomputed(he, block) ) {
              LOG << "Pin count in part of hyperedge" << he << "in block" << block
                  << "is" << pinCountInPart(he, block) << ", but should be"
                  << pinCountInPartRecomputed(he, block);
              return false;
            }
          }
          return true;
        }());
      }
    });
  }

  // ####################### Partition Information #######################

  // ! Block that vertex u belongs to
  PartitionID partID(const HypernodeID u) const {
    ASSERT(u < initialNumNodes(), "Hypernode" << u << "does not exist");
    return _part_ids[u];
  }

  void extractPartIDs(Array<PartitionID>& part_ids) {
    std::swap(_part_ids, part_ids);
  }

  void setOnlyNodePart(const HypernodeID u, PartitionID p) {
    ASSERT(p != kInvalidPartition && p < _k);
    ASSERT(_part_ids[u] == kInvalidPartition);
    _part_ids[u] = p;
  }

  void setNodePart(const HypernodeID u, PartitionID p) {
    setOnlyNodePart(u, p);
    _part_weights[p].fetch_add(nodeWeight(u), std::memory_order_relaxed);
    for (HyperedgeID he : incidentEdges(u)) {
      incrementPinCountInPartWithoutGainUpdate(he, p);
    }
  }

  // ! Changes the block id of vertex u from block 'from' to block 'to'
  // ! Returns true, if move of vertex u to corresponding block succeeds.
  template<typename SuccessFunc, typename DeltaFunc>
  bool changeNodePart(const HypernodeID u,
                      PartitionID from,
                      PartitionID to,
                      HypernodeWeight max_weight_to,
                      SuccessFunc&& report_success,
                      DeltaFunc&& delta_func) {
    ASSERT(partID(u) == from);
    ASSERT(from != to);
    const HypernodeWeight wu = nodeWeight(u);
    const HypernodeWeight to_weight_after = _part_weights[to].add_fetch(wu, std::memory_order_relaxed);
    if (to_weight_after <= max_weight_to) {
      _part_ids[u] = to;
      _part_weights[from].fetch_sub(wu, std::memory_order_relaxed);
      report_success();
      for ( const HyperedgeID he : incidentEdges(u) ) {
        updatePinCountOfHyperedge(he, from, to, delta_func);
      }
      return true;
    } else {
      _part_weights[to].fetch_sub(wu, std::memory_order_relaxed);
      return false;
    }
  }

  // curry
  bool changeNodePart(const HypernodeID u,
                      PartitionID from,
                      PartitionID to,
                      const DeltaFunction& delta_func = NOOP_FUNC) {
    return changeNodePart(u, from, to,
      std::numeric_limits<HypernodeWeight>::max(), []{}, delta_func);
  }

  // Make sure not to call phg.gainCacheUpdate(..) in delta_func for changeNodePartWithGainCacheUpdate
  template<typename SuccessFunc, typename DeltaFunc>
  bool changeNodePartWithGainCacheUpdate(const HypernodeID u,
                                         PartitionID from,
                                         PartitionID to,
                                         HypernodeWeight max_weight_to,
                                         SuccessFunc&& report_success,
                                         DeltaFunc&& delta_func) {
    //ASSERT(_is_gain_cache_initialized, "Gain cache is not initialized");

    auto my_delta_func = [&](const HyperedgeID he, const HyperedgeWeight edge_weight, const HypernodeID edge_size,
            const HypernodeID pin_count_in_from_part_after, const HypernodeID pin_count_in_to_part_after) {
      delta_func(he, edge_weight, edge_size, pin_count_in_from_part_after, pin_count_in_to_part_after);
      gainCacheUpdate(he, edge_weight, from, pin_count_in_from_part_after, to, pin_count_in_to_part_after);
    };
    return changeNodePart(u, from, to, max_weight_to, report_success, my_delta_func);

  }

  bool changeNodePartWithGainCacheUpdate(const HypernodeID u, PartitionID from, PartitionID to) {
    return changeNodePartWithGainCacheUpdate(u, from, to, std::numeric_limits<HypernodeWeight>::max(), [] { },
                                             NoOpDeltaFunc());
  }

  // ! Weight of a block
  HypernodeWeight partWeight(const PartitionID p) const {
    ASSERT(p != kInvalidPartition && p < _k);
    return _part_weights[p].load(std::memory_order_relaxed);
  }

  // ! Returns, whether hypernode u is adjacent to a least one cut hyperedge.
  bool isBorderNode(const HypernodeID u) const {
    if ( nodeDegree(u) <= HIGH_DEGREE_THRESHOLD ) {
      for ( const HyperedgeID& he : incidentEdges(u) ) {
        if ( connectivity(he) > 1 ) {
          return true;
        }
      }
      return false;
    } else {
      // TODO maybe we should allow these in label propagation? definitely not in FM
      // In case u is a high degree vertex, we omit the border node check and
      // and return false. Assumption is that it is very unlikely that such a
      // vertex can change its block.
      return false;
    }
  }

  HypernodeID numIncidentCutHyperedges(const HypernodeID u) const {
    HypernodeID num_incident_cut_hyperedges = 0;
    for ( const HyperedgeID& he : incidentEdges(u) ) {
      if ( connectivity(he) > 1 ) {
        ++num_incident_cut_hyperedges;
      }
    }
    return num_incident_cut_hyperedges;
  }

  // ! Number of blocks which pins of hyperedge e belongs to
  PartitionID connectivity(const HyperedgeID e) const {
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    return _connectivity_set.connectivity(e);
  }

  // ! Returns the number pins of hyperedge e that are part of block id
  HypernodeID pinCountInPart(const HyperedgeID e, const PartitionID p) const {
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(p != kInvalidPartition && p < _k);
    return _pins_in_part.pinCountInPart(e, p);
  }

  /**
   * In the following, we define functions to derive the benefit and penalty term
   * decribed in our publications. However, we swapped the naming of both (the benefit term
   * is now the penalty term and vice versa) due to refactoring of the gain cache.
   */

  // ! The move from penalty term stores weight of all incident edges of u for which
  // ! we cannot reduce their connecitivity when u is moved out of its block.
  // ! More formally, p(u) := w({ e \in I(u) | pin_count(e, partID(u)) > 1 })
  HyperedgeWeight moveFromPenalty(const HypernodeID u) const {
    //ASSERT(_is_gain_cache_initialized, "Gain cache is not initialized");
    return _gain_cache[penalty_index(u)].load(std::memory_order_relaxed);
  }

  // ! The move to benefit term stores the weight of all incident edges of u
  // ! for which we do not increase the connecitivity when we move u to block p.
  // ! More formally, b(u, p) := w({ e \in I(u) | pin_count(e, p) >= 1 })
  HyperedgeWeight moveToBenefit(const HypernodeID u, PartitionID p) const {
    //ASSERT(_is_gain_cache_initialized, "Gain cache is not initialized");
    return _gain_cache[benefit_index(u, p)].load(std::memory_order_relaxed);
  }

  // ! The gain of moving a node u from its current block to a target block to can
  // ! be expressed as g(u, to) = b(u, to) - p(u), which is the weight of
  // ! all incident edges of u for which we do not increase their connectivity if moved
  // ! to block to minus the weight of all incident edges of u for which we cannot reduce
  // ! their connectivity if u is moved out of its current block.
  HyperedgeWeight km1Gain(const HypernodeID u, PartitionID from, PartitionID to) const {
    unused(from);
    //ASSERT(_is_gain_cache_initialized, "Gain cache is not initialized");
    ASSERT(from == partID(u), "While gain computation works for from != partID(u), such a query makes no sense");
    ASSERT(from != to, "The gain computation doesn't work for from = to");
    return moveToBenefit(u, to) - moveFromPenalty(u);
  }

  void allocateGainTableIfNecessary() {
    if (_gain_cache.size() == 0) {
      _gain_cache.resize(
              "Refinement", "gain_cache", _top_level_num_nodes * size_t(_k + 1), true);
    }
  }


  // ! The move from penalty term stores weight of all incident edges of u for which
  // ! we cannot reduce their connecitivity when u is moved out of its block.
  // ! More formally, p(u) := w({ e \in I(u) | pin_count(e, partID(u)) > 1 })
  HyperedgeWeight moveFromPenaltyRecomputed(const HypernodeID u) const {
    const PartitionID p = partID(u);
    HyperedgeWeight b = 0;
    for (HyperedgeID e : incidentEdges(u)) {
      if (pinCountInPart(e, p) > 1) {
        b += edgeWeight(e);
      }
    }
    return b;
  }

  void recomputeMoveFromPenalty(const HypernodeID u) {
    _gain_cache[penalty_index(u)].store(moveFromPenaltyRecomputed(u), std::memory_order_relaxed);
  }

  // ! Only for testing
  // ! The move to benefit term stores the weight of all incident edges of u
  // ! for which we do not increase the connecitivity when we move u to block p.
  // ! More formally, b(u, p) := w({ e \in I(u) | pin_count(e, p) >= 1 })
  HyperedgeWeight moveToBenefitRecomputed(const HypernodeID u, PartitionID p) const {
    HyperedgeWeight w = 0;
    for (HyperedgeID e : incidentEdges(u)) {
      if (pinCountInPart(e, p) >= 1) {
        w += edgeWeight(e);
      }
    }
    return w;
  }

  // ! Initializes the partition of the hypergraph, if block ids are assigned with
  // ! setOnlyNodePart(...). In that case, block weights and pin counts in part for
  // ! each hyperedge must be initialized explicitly here.
  void initializePartition() {
    tbb::parallel_invoke(
            [&] { initializeBlockWeights(); },
            [&] { initializePinCountInPart(); }
    );
  }

  bool isGainCacheInitialized() const {
    return _is_gain_cache_initialized;
  }

  void initializeGainCacheEntry(const HypernodeID u, vec<Gain>& benefit_aggregator) {
    PartitionID pu = partID(u);
    Gain penalty = 0;
    for (const HyperedgeID& e : incidentEdges(u)) {
      HyperedgeWeight ew = edgeWeight(e);
      if (pinCountInPart(e, pu) > 1) {
        penalty += ew;
      }
      for (const PartitionID& i : connectivitySet(e)) {
        benefit_aggregator[i] += ew;
      }
    }

    _gain_cache[penalty_index(u)].store(penalty, std::memory_order_relaxed);
    for (PartitionID i = 0; i < _k; ++i) {
      _gain_cache[benefit_index(u, i)].store(benefit_aggregator[i], std::memory_order_relaxed);
      benefit_aggregator[i] = 0;
    }
  }

  // ! Initialize gain cache
  // ! NOTE: Requires that pin counts are already initialized and reflect the
  // ! current state of the partition
  void initializeGainCache() {
    allocateGainTableIfNecessary();

    // check whether part has been initialized
    ASSERT(std::none_of(nodes().begin(), nodes().end(),
                            [&](HypernodeID u) { return partID(u) == kInvalidPartition || partID(u) > k(); }) );

    auto aggregate_contribution_of_he_for_vertex =
      [&](const PartitionID block_of_u,
          const HyperedgeID he,
          HyperedgeWeight& l_move_from_penalty,
          vec<HyperedgeWeight>& l_move_to_benefit) {
      HyperedgeWeight edge_weight = edgeWeight(he);
      if (pinCountInPart(he, block_of_u) > 1) {
        l_move_from_penalty += edge_weight;
      }

      for (const PartitionID block : connectivitySet(he)) {
        l_move_to_benefit[block] += edge_weight;
      }
    };



    // Gain calculation consist of two stages
    //  1. Compute gain of all low degree vertices sequential (iterating over all vertices in parallel)
    //  2. Compute gain of all high degree vertices parallel (iterating over all high degree vertices sequentially)
    tbb::enumerable_thread_specific< vec<HyperedgeWeight> > ets_mtb(_k, 0);
    std::mutex high_degree_vertex_mutex;
    parallel::scalable_vector<HypernodeID> high_degree_vertices;

    // Compute gain of all low degree vertices sequential (iterating over all vertices in parallel)
    tbb::parallel_for(tbb::blocked_range<HypernodeID>(HypernodeID(0), initialNumNodes()),
      [&](tbb::blocked_range<HypernodeID>& r) {
        vec<HyperedgeWeight>& l_move_to_benefit = ets_mtb.local();
        for (HypernodeID u = r.begin(); u < r.end(); ++u) {
          if ( nodeIsEnabled(u)) {
            if ( nodeDegree(u) <= HIGH_DEGREE_THRESHOLD) {
              const PartitionID from = partID(u);
              HyperedgeWeight l_move_from_penalty = 0;
              for (HyperedgeID he : incidentEdges(u)) {
                aggregate_contribution_of_he_for_vertex(from, he,
                  l_move_from_penalty, l_move_to_benefit);
              }

              _gain_cache[penalty_index(u)].store(l_move_from_penalty, std::memory_order_relaxed);
              for (PartitionID p = 0; p < _k; ++p) {
                _gain_cache[benefit_index(u,p)].store(l_move_to_benefit[p], std::memory_order_relaxed);
                l_move_to_benefit[p] = 0;
              }
            } else {
              // Collect high degree vertex for subsequent parallel gain computation
              std::lock_guard<std::mutex> lock(high_degree_vertex_mutex);
              high_degree_vertices.push_back(u);
            }
          }
        }
      });

    // Compute gain of all high degree vertices parallel (iterating over all high degree vertices sequentially)
    for ( const HypernodeID& u : high_degree_vertices ) {
      tbb::enumerable_thread_specific<HyperedgeWeight> ets_mfp(0);
      const PartitionID from = partID(u);
      const HypernodeID degree_of_u = _hg->nodeDegree(u);
      tbb::parallel_for(tbb::blocked_range<HypernodeID>(ID(0), degree_of_u),
        [&](tbb::blocked_range<HypernodeID>& r) {
        vec<HyperedgeWeight>& l_move_to_benefit = ets_mtb.local();
        HyperedgeWeight& l_move_from_penalty = ets_mfp.local();
        size_t current_pos = r.begin();
        for ( const HyperedgeID& he : _hg->incident_nets_of(u, r.begin()) ) {
          aggregate_contribution_of_he_for_vertex(from, he,
            l_move_from_penalty, l_move_to_benefit);
          ++current_pos;
          if ( current_pos == r.end() ) {
            break;
          }
        }
      });

      // Aggregate thread locals to compute overall gain of the high degree vertex
      const HyperedgeWeight penalty_term = ets_mfp.combine(std::plus<HyperedgeWeight>());
      _gain_cache[penalty_index(u)].store(penalty_term, std::memory_order_relaxed);
      for (PartitionID p = 0; p < _k; ++p) {
        HyperedgeWeight move_to_benefit = 0;
        for ( auto& l_move_to_benefit : ets_mtb ) {
          move_to_benefit += l_move_to_benefit[p];
          l_move_to_benefit[p] = 0;
        }
        _gain_cache[benefit_index(u, p)].store(move_to_benefit, std::memory_order_relaxed);
      }
    }

    _is_gain_cache_initialized = true;
  }

  // ! Reset partition (not thread-safe)
  void resetPartition() {
    _part_ids.assign(_part_ids.size(), kInvalidPartition, false);
    for (auto& x : _part_weights) x.store(0, std::memory_order_relaxed);

    // Reset pin count in part and connectivity set
    for ( const HyperedgeID& he : edges() ) {
      for ( const PartitionID& block : connectivitySet(he) ) {
        _pins_in_part.setPinCountInPart(he, block, 0);
      }
      _connectivity_set.clear(he);
    }
  }

  // ! Should be called e.g. after a rollback (see PartitionedGraph).
  void resetMoveState() {
    // Nothing to do here
  }

  // ! Only for testing
  void recomputePartWeights() {
    for (PartitionID p = 0; p < _k; ++p) {
      _part_weights[p].store(0);
    }

    for (HypernodeID u : nodes()) {
      _part_weights[ partID(u) ] += nodeWeight(u);
    }
  }

  // ! Only for testing
  bool checkTrackedPartitionInformation() {
    bool success = true;

    for (HyperedgeID e : edges()) {
      PartitionID expected_connectivity = 0;
      for (PartitionID i = 0; i < k(); ++i) {
        const HypernodeID actual_pin_count_in_part = pinCountInPart(e, i);
        if ( actual_pin_count_in_part != pinCountInPartRecomputed(e, i) ) {
          LOG << "Pin count of hyperedge" << e << "in block" << i << "=>" <<
              "Expected:" << V(pinCountInPartRecomputed(e, i)) << "," <<
              "Actual:" <<  V(pinCountInPart(e, i));
          success = false;
        }
        expected_connectivity += (actual_pin_count_in_part > 0);
      }
      if ( expected_connectivity != connectivity(e) ) {
        LOG << "Connectivity of hyperedge" << e << "=>" <<
            "Expected:" << V(expected_connectivity)  << "," <<
            "Actual:" << V(connectivity(e));
        success = false;
      }
    }

    if ( _is_gain_cache_initialized ) {
      for (HypernodeID u : nodes()) {
        if ( moveFromPenalty(u) != moveFromPenaltyRecomputed(u) ) {
          LOG << "Move from benefit of hypernode" << u << "=>" <<
              "Expected:" << V(moveFromPenaltyRecomputed(u)) << ", " <<
              "Actual:" <<  V(moveFromPenalty(u));
          for ( const HyperedgeID& e : incidentEdges(u) ) {
            LOG << V(u) << V(partID(u)) << V(e) << V(edgeSize(e)) << V(edgeWeight(e)) << V(pinCountInPart(e, partID(u)));
          }
          success = false;
        }

        for (PartitionID i = 0; i < k(); ++i) {
          if (partID(u) != i) {
            if ( moveToBenefit(u, i) != moveToBenefitRecomputed(u, i) ) {
              LOG << "Move to penalty of hypernode" << u << "in block" << i << "=>" <<
                  "Expected:" << V(moveToBenefitRecomputed(u, i)) << ", " <<
                  "Actual:" <<  V(moveToBenefit(u, i));
              success = false;
            }
          }
        }
      }
    }
    return success;
  }

  // ####################### Memory Consumption #######################

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    utils::MemoryTreeNode* hypergraph_node = parent->addChild("Hypergraph");
    _hg->memoryConsumption(hypergraph_node);
    utils::MemoryTreeNode* connectivity_set_node = parent->addChild("Connectivity Sets");
    _connectivity_set.memoryConsumption(connectivity_set_node);

    parent->addChild("Part Weights", sizeof(CAtomic<HypernodeWeight>) * _k);
    parent->addChild("Part IDs", sizeof(PartitionID) * _hg->initialNumNodes());
    parent->addChild("Pin Count In Part", _pins_in_part.size_in_bytes());
    parent->addChild("Gain Cache", sizeof(HyperedgeWeight) * _gain_cache.size());
    parent->addChild("HE Ownership", sizeof(SpinLock) * _hg->initialNumNodes());
  }

  // ####################### Extract Block #######################

  // ! Extracts a block of a partition as separate hypergraph.
  // ! It also returns a vertex-mapping from the original hypergraph to the sub-hypergraph.
  // ! If cut_net_splitting is activated, hyperedges that span more than one block (cut nets) are split, which is used for the connectivity metric.
  // ! Otherwise cut nets are discarded (cut metric).
  std::pair<Hypergraph, parallel::scalable_vector<HypernodeID> > extract(
          PartitionID block,
          bool cut_net_splitting,
          bool stable_construction_of_incident_edges) {
    ASSERT(block != kInvalidPartition && block < _k);

    // Compactify vertex ids
    parallel::scalable_vector<HypernodeID> hn_mapping(_hg->initialNumNodes(), kInvalidHypernode);
    parallel::scalable_vector<HyperedgeID> he_mapping(_hg->initialNumEdges(), kInvalidHyperedge);
    HypernodeID num_hypernodes = 0;
    HypernodeID num_hyperedges = 0;
    tbb::parallel_invoke([&] {
      for ( const HypernodeID& hn : nodes() ) {
        if ( partID(hn) == block ) {
          hn_mapping[hn] = num_hypernodes++;
        }
      }
    }, [&] {
      for ( const HyperedgeID& he : edges() ) {
        if ( pinCountInPart(he, block) > 1 &&
             (cut_net_splitting || connectivity(he) == 1) ) {
          he_mapping[he] = num_hyperedges++;
        }
      }
    });

    // Extract plain hypergraph data for corresponding block
    using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;
    HyperedgeVector edge_vector;
    parallel::scalable_vector<HyperedgeWeight> hyperedge_weight;
    parallel::scalable_vector<HypernodeWeight> hypernode_weight;
    tbb::parallel_invoke([&] {
      edge_vector.resize(num_hyperedges);
      hyperedge_weight.resize(num_hyperedges);
      doParallelForAllEdges([&](const HyperedgeID he) {
        if ( pinCountInPart(he, block) > 1 &&
             (cut_net_splitting || connectivity(he) == 1) ) {
          ASSERT(he_mapping[he] < num_hyperedges);
          hyperedge_weight[he_mapping[he]] = edgeWeight(he);
          for ( const HypernodeID& pin : pins(he) ) {
            if ( partID(pin) == block ) {
              edge_vector[he_mapping[he]].push_back(hn_mapping[pin]);
            }
          }
        }
      });
    }, [&] {
      hypernode_weight.resize(num_hypernodes);
      doParallelForAllNodes([&](const HypernodeID hn) {
        if ( partID(hn) == block ) {
          hypernode_weight[hn_mapping[hn]] = nodeWeight(hn);
        }
      });
    });

    // Construct hypergraph
    Hypergraph extracted_hypergraph = HypergraphFactory::construct(num_hypernodes, num_hyperedges,
            edge_vector, hyperedge_weight.data(), hypernode_weight.data(), stable_construction_of_incident_edges);

    // Set community ids
    doParallelForAllNodes([&](const HypernodeID& hn) {
      if ( partID(hn) == block ) {
        const HypernodeID extracted_hn = hn_mapping[hn];
        extracted_hypergraph.setCommunityID(extracted_hn, _hg->communityID(hn));
      }
    });
    return std::make_pair(std::move(extracted_hypergraph), std::move(hn_mapping));
  }

  void freeInternalData() {
    if ( _k > 0 ) {
      tbb::parallel_invoke( [&] {
        parallel::parallel_free(_part_ids, _pin_count_update_ownership);
      }, [&] {
        parallel::free(_pins_in_part.data());
      }, [&] {
        _connectivity_set.freeInternalData();
      } );
    }
    _k = 0;
  }


  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void gainCacheUpdate(const HyperedgeID he, const HyperedgeWeight we,
                       const PartitionID from, const HypernodeID pin_count_in_from_part_after,
                       const PartitionID to, const HypernodeID pin_count_in_to_part_after) {
    if (pin_count_in_from_part_after == 1) {
      for (const HypernodeID& u : pins(he)) {
        nodeGainAssertions(u, from);
        if (partID(u) == from) {
          _gain_cache[penalty_index(u)].fetch_sub(we, std::memory_order_relaxed);
        }
      }
    } else if (pin_count_in_from_part_after == 0) {
      for (const HypernodeID& u : pins(he)) {
        nodeGainAssertions(u, from);
        _gain_cache[benefit_index(u, from)].fetch_sub(we, std::memory_order_relaxed);
      }
    }

    if (pin_count_in_to_part_after == 1) {
      for (const HypernodeID& u : pins(he)) {
        nodeGainAssertions(u, to);
        _gain_cache[benefit_index(u, to)].fetch_add(we, std::memory_order_relaxed);
      }
    } else if (pin_count_in_to_part_after == 2) {
      for (const HypernodeID& u : pins(he)) {
        nodeGainAssertions(u, to);
        if (partID(u) == to) {
          _gain_cache[penalty_index(u)].fetch_add(we, std::memory_order_relaxed);
        }
      }
    }
  }

 private:

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  size_t penalty_index(const HypernodeID u) const {
    return size_t(u) * ( _k + 1 );
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  size_t benefit_index(const HypernodeID u, const PartitionID p) const {
    return size_t(u) * ( _k + 1 )  + p + 1;
  }

  void applyPartWeightUpdates(vec<HypernodeWeight>& part_weight_deltas) {
    for (PartitionID p = 0; p < _k; ++p) {
      _part_weights[p].fetch_add(part_weight_deltas[p], std::memory_order_relaxed);
    }
  }

  void initializeBlockWeights() {
    auto accumulate = [&](tbb::blocked_range<HypernodeID>& r) {
      vec<HypernodeWeight> pws(_k, 0);  // this is not enumerable_thread_specific because of the static partitioner
      for (HypernodeID u = r.begin(); u < r.end(); ++u) {
        if ( nodeIsEnabled(u) ) {
          const PartitionID pu = partID( u );
          const HypernodeWeight wu = nodeWeight( u );
          pws[pu] += wu;
        }
      }
      applyPartWeightUpdates(pws);
    };

    tbb::parallel_for(tbb::blocked_range<HypernodeID>(HypernodeID(0), initialNumNodes()),
                      accumulate,
                      tbb::static_partitioner()
    );
  }

  void initializePinCountInPart() {
    tls_enumerable_thread_specific< vec<HypernodeID> > ets_pin_count_in_part(_k, 0);

    auto assign = [&](tbb::blocked_range<HyperedgeID>& r) {
      vec<HypernodeID>& pin_counts = ets_pin_count_in_part.local();
      for (HyperedgeID he = r.begin(); he < r.end(); ++he) {
        if ( edgeIsEnabled(he) ) {
          for (const HypernodeID& pin : pins(he)) {
            ++pin_counts[partID(pin)];
          }

          for (PartitionID p = 0; p < _k; ++p) {
            ASSERT(pinCountInPart(he, p) == 0);
            if (pin_counts[p] > 0) {
              _connectivity_set.add(he, p);
              _pins_in_part.setPinCountInPart(he, p, pin_counts[p]);
            }
            pin_counts[p] = 0;
          }
        }
      }
    };

    tbb::parallel_for(tbb::blocked_range<HyperedgeID>(HyperedgeID(0), initialNumEdges()), assign);
  }

  HypernodeID pinCountInPartRecomputed(const HyperedgeID e, PartitionID p) const {
    HypernodeID pcip = 0;
    for (HypernodeID u : pins(e)) {
      if (partID(u) == p) {
        pcip++;
      }
    }
    return pcip;
  }

  void nodeGainAssertions(const HypernodeID u, const PartitionID p) const {
    unused(u);
    unused(p);
    ASSERT(u < initialNumNodes(), "Hypernode" << u << "does not exist");
    ASSERT(nodeIsEnabled(u), "Hypernode" << u << "is disabled");
    ASSERT(p != kInvalidPartition && p < _k);
    ASSERT(benefit_index(u, p) < _gain_cache.size());
  }

  // ! Updates pin count in part using a spinlock.
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void updatePinCountOfHyperedge(const HyperedgeID he,
                                                                    const PartitionID from,
                                                                    const PartitionID to,
                                                                    const DeltaFunction& delta_func) {
    ASSERT(he < _pin_count_update_ownership.size());
    _pin_count_update_ownership[he].lock();
    const HypernodeID pin_count_in_from_part_after = decrementPinCountInPartWithoutGainUpdate(he, from);
    const HypernodeID pin_count_in_to_part_after = incrementPinCountInPartWithoutGainUpdate(he, to);
    _pin_count_update_ownership[he].unlock();
    delta_func(he, edgeWeight(he), edgeSize(he), pin_count_in_from_part_after, pin_count_in_to_part_after);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HypernodeID decrementPinCountInPartWithoutGainUpdate(const HyperedgeID e, const PartitionID p) {
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(p != kInvalidPartition && p < _k);
    const HypernodeID pin_count_after = _pins_in_part.decrementPinCountInPart(e, p);
    if ( pin_count_after == 0 ) {
      _connectivity_set.remove(e, p);
    }
    return pin_count_after;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HypernodeID incrementPinCountInPartWithoutGainUpdate(const HyperedgeID e, const PartitionID p) {
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(p != kInvalidPartition && p < _k);
    const HypernodeID pin_count_after = _pins_in_part.incrementPinCountInPart(e, p);
    if ( pin_count_after == 1 ) {
      _connectivity_set.add(e, p);
    }
    return pin_count_after;
  }


  // ! Indicate wheater gain cache is initialized
  bool _is_gain_cache_initialized;

  size_t _top_level_num_nodes = 0;

  // ! Number of blocks
  PartitionID _k = 0;

  // ! Hypergraph object around which this partitioned hypergraph is wrapped
  Hypergraph* _hg = nullptr;

  // ! Weight and information for all blocks.
  vec< CAtomic<HypernodeWeight> > _part_weights;

  // ! Current block IDs of the vertices
  Array< PartitionID > _part_ids;

  // ! For each hyperedge and each block, _pins_in_part stores the
  // ! number of pins in that block
  PinCountInPart _pins_in_part;

  // ! For each hyperedge, _connectivity_set stores the set of blocks that the hyperedge spans
  ConnectivitySets _connectivity_set;

  // ! The gain of moving a node u to from its current block V_i to a target block V_j
  // ! can be expressed as follows for the connectivity metric
  // ! g(u, V_j) := w({ e \in I(u) | pin_count(e, V_j) >= 1 }) - w({ e \in I(u) | pin_count(e, V_i) > 1 })
  // !            = b(u, V_j) - p(u)
  // ! We call b(u, V_j) the benefit term and p(u) the penalty term. Our gain cache stores and maintains these
  // ! entries for each node and block. Thus, the gain cache stores k + 1 entries per node.
  Array< CAtomic<HyperedgeWeight> > _gain_cache;

  // ! In order to update the pin count of a hyperedge thread-safe, a thread must acquire
  // ! the ownership of a hyperedge via a CAS operation.
  Array<SpinLock> _pin_count_update_ownership;
};

} // namespace ds
} // namespace mt_kahypar

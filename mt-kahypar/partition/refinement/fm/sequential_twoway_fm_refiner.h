/*******************************************************************************
 * This file is part of MT-KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "kahypar/datastructure/kway_priority_queue.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"

namespace mt_kahypar {

/**
 * Implements a classical sequential 2-way FM which is similiar to the one implemented in KaHyPar.
 * Main difference is that we do not use a gain cache, since we do not want use the 2-way fm refiner
 * in a multilevel context. It is used after each bisection during initial partitioning to refine
 * a given bipartition.
 */
class SequentialTwoWayFmRefiner {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  using KWayRefinementPQ = kahypar::ds::KWayPriorityQueue<HypernodeID, Gain, std::numeric_limits<Gain> >;

  /**
   * A hyperedge can be in three states during FM refinement: FREE, LOOSE and LOCKED.
   * Initial all hyperedges are FREE. Once we move a vertex incident to the hyperedge,
   * the hyperedge becomes LOOSE. If we move an other vertex in the opposite direction
   * the hyperedge becomes LOCKED. LOCKED hyperedges have the property that they can not
   * be removed from cut and we can therefore skip delta gain updates.
   */
  enum HEState {
    FREE = std::numeric_limits<PartitionID>::max() - 1,
    LOCKED = std::numeric_limits<PartitionID>::max(),
  };

  /**
   * INACTIVE = Initial state of a vertex
   * ACTIVE = Vertex is a border node and inserted into the PQ
   * MOVED = Vertex was moved during local search
   */
  enum class VertexState {
    INACTIVE,
    ACTIVE,
    MOVED
  };

  /**
   * Our partitioned hypergraph data structures does not track to how many cut hyperedges
   * a vertex is incident to. This helper class tracks the number of incident cut hyperedges
   * and gathers all nodes that became border or internal nodes during the last move on
   * the hypergraph, which is required by our 2-way FM implementation, since only border vertices
   * are eligble for moving.
   */
  class BorderVertexTracker {

   public:
    explicit BorderVertexTracker(const HypernodeID& num_hypernodes) :
      _num_hypernodes(num_hypernodes),
      _num_incident_cut_hes(num_hypernodes, 0),
      _hns_to_activate(),
      _hns_to_remove_from_pq() { }

    void initialize(const PartitionedHypergraph& phg) {
      reset();
      for ( const HypernodeID& hn : phg.nodes() ) {
        ASSERT(hn <  _num_hypernodes);
        for ( const HyperedgeID& he : phg.incidentEdges(hn) ) {
          if ( phg.connectivity(he) > 1 ) {
            ++_num_incident_cut_hes[hn];
          }
        }
      }
    }

    bool isBorderNode(const HypernodeID hn) const {
      ASSERT(hn <  _num_hypernodes);
      return _num_incident_cut_hes[hn] > 0;
    }

    void becameCutHyperedge(const PartitionedHypergraph& phg,
                            const HyperedgeID he,
                            const parallel::scalable_vector<VertexState>& vertex_state) {
      ASSERT(phg.connectivity(he) > 1);
      for ( const HypernodeID& pin : phg.pins(he) ) {
        ASSERT(pin <  _num_hypernodes);
        ASSERT(_num_incident_cut_hes[pin] <= phg.nodeDegree(pin));
        ++_num_incident_cut_hes[pin];
        if ( _num_incident_cut_hes[pin] == 1 && vertex_state[pin] == VertexState::INACTIVE ) {
          _hns_to_activate.push_back(pin);
        }
      }
    }

    void becameNonCutHyperedge(const PartitionedHypergraph& phg,
                               const HyperedgeID he,
                               const parallel::scalable_vector<VertexState>& vertex_state) {
      ASSERT(phg.connectivity(he) == 1);
      for ( const HypernodeID& pin : phg.pins(he) ) {
        ASSERT(pin <  _num_hypernodes);
        ASSERT(_num_incident_cut_hes[pin] > 0);
        --_num_incident_cut_hes[pin];
        // Note, it is possible that we insert border vertices here, since an other hyperedge
        // can become cut after the update. However, we handle this later by an explicit check
        // if the vertex is still an internal vertex (see doForAllVerticesThatBecameInternalVertices(...)).
        if ( _num_incident_cut_hes[pin] == 0 && vertex_state[pin] == VertexState::ACTIVE ) {
          _hns_to_remove_from_pq.push_back(pin);
        }
      }
    }

    // ! Iterates over all vertices that became border vertices after the last move
    template<typename F>
    void doForAllVerticesThatBecameBorderVertices(const F& f) {
      for ( const HypernodeID& hn : _hns_to_activate ) {
        f(hn);
      }
      _hns_to_activate.clear();
    }

    // ! Iterates over all vertices that became internal vertices after the last move
    template<typename F>
    void doForAllVerticesThatBecameInternalVertices(const F& f) {
      for ( const HypernodeID& hn : _hns_to_remove_from_pq ) {
        // Explicit border vertex check, because vector can contain fales positives
        // (see becameNonCutHyperedge(...))
        if ( !isBorderNode(hn) ) {
          f(hn);
        }
      }
      _hns_to_remove_from_pq.clear();
    }

   private:
    void reset() {
      for ( HypernodeID hn = 0; hn < _num_hypernodes; ++hn ) {
        _num_incident_cut_hes[hn] = 0;
      }
      _hns_to_activate.clear();
      _hns_to_remove_from_pq.clear();
    }

    const HypernodeID _num_hypernodes;
    parallel::scalable_vector<HyperedgeID> _num_incident_cut_hes;
    parallel::scalable_vector<HypernodeID> _hns_to_activate;
    parallel::scalable_vector<HypernodeID> _hns_to_remove_from_pq;
  };

 public:
  SequentialTwoWayFmRefiner(PartitionedHypergraph& phg,
                            const Context& context) :
    _phg(phg),
    _context(context),
    _nodes(),
    _pq(context.partition.k),
    _border_vertices(phg.initialNumNodes()),
    _vertex_state(phg.initialNumNodes(), VertexState::INACTIVE),
    _he_state(phg.initialNumEdges(), HEState::FREE) {
    ASSERT(_context.partition.k == 2);
    _pq.initialize(_phg.initialNumNodes());
    for ( const HypernodeID& hn : phg.nodes() ) {
      _nodes.push_back(hn);
    }
  }

  bool refine(kahypar::Metrics& best_metrics);

 private:

  void activate(const HypernodeID hn);

  /**
   * Performs delta gain update on all non locked hyperedges and
   * state transition of hyperedges.
   */
  void updateNeighbors(const HypernodeID hn,
                       const PartitionID from,
                       const PartitionID to);

  // ! Delta-Gain Update as decribed in [ParMar06].
  void deltaGainUpdate(const HyperedgeID he,
                       const PartitionID from,
                       const PartitionID to);

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void updatePin(const HypernodeID pin, const Gain delta);

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void updatePQState(const PartitionID from,
                                                        const PartitionID to);

  Gain computeGain(const HypernodeID hn, const PartitionID from, const PartitionID to);

  void rollback(const parallel::scalable_vector<HypernodeID>& performed_moves,
                const size_t min_cut_idx);

  bool verifyPQState() const;

  PartitionedHypergraph& _phg;
  const Context& _context;

  parallel::scalable_vector<HypernodeID> _nodes;
  KWayRefinementPQ _pq;
  BorderVertexTracker _border_vertices;
  parallel::scalable_vector<VertexState> _vertex_state;
  parallel::scalable_vector<PartitionID> _he_state;
};

} // namespace mt_kahypar

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
#include "kahypar/partition/metrics.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/fm/stop_rule.h"


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
   * TODO(heuer): We should integrate this again into the partitioned hypergraph, since we should
   * also change our parallel k-way implementation to only move border vertices.
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
    _pq(context.partition.k),
    _border_vertices(phg.initialNumNodes()),
    _vertex_state(phg.initialNumNodes(), VertexState::INACTIVE),
    _he_state(phg.initialNumEdges(), HEState::FREE) {
    ASSERT(_context.partition.k == 2);
    _pq.initialize(_phg.initialNumNodes());
  }

  bool refine(kahypar::Metrics& best_metrics) {

    // Activate all border nodes
    _pq.clear();
    _border_vertices.initialize(_phg);
    for ( const HypernodeID& hn : _phg.nodes() ) {
      _vertex_state[hn] = VertexState::INACTIVE;
      activate(hn);
    }
    for ( const HyperedgeID& he : _phg.edges() ) {
      _he_state[he] = HEState::FREE;
    }

    auto border_vertex_update = [&](const HyperedgeID he,
                                    const HyperedgeWeight,
                                    const HypernodeID edge_size,
                                    const HypernodeID pin_count_in_from_part_after,
                                    const HypernodeID pin_count_in_to_part_after) {
                             if ( edge_size > 1 ) {
                               if ( pin_count_in_from_part_after == 0 ) {
                                 _border_vertices.becameNonCutHyperedge(_phg, he, _vertex_state);
                               } else if ( pin_count_in_to_part_after == 1 ) {
                                 _border_vertices.becameCutHyperedge(_phg, he, _vertex_state);
                               }
                             }
                           };

    parallel::scalable_vector<HypernodeID> performed_moves;
    HyperedgeWeight current_cut = best_metrics.cut;
    double current_imbalance = best_metrics.imbalance;
    size_t min_cut_idx = 0;
    StopRule stopping_rule(_phg.initialNumNodes());
    while ( !_pq.empty() && !stopping_rule.searchShouldStop() ) {
      ASSERT(_pq.isEnabled(0) || _pq.isEnabled(1));
      HEAVY_REFINEMENT_ASSERT(verifyPQState(), "PQ corrupted!");

      // Retrieve max gain move from PQ
      Gain gain = invalidGain;
      HypernodeID hn = kInvalidHypernode;
      PartitionID to = kInvalidPartition;
      _pq.deleteMax(hn, gain, to);

      ASSERT(hn != kInvalidHypernode);
      ASSERT(_border_vertices.isBorderNode(hn));
      ASSERT(_phg.partID(hn) == 1 - to);
      HEAVY_REFINEMENT_ASSERT(gain == computeGain(hn, _phg.partID(hn), to));

      // Perform vertex move
      PartitionID from = _phg.partID(hn);
      _vertex_state[hn] = VertexState::MOVED;
      if ( _phg.changeNodePart(hn, from, to,
            _context.partition.max_part_weights[to], border_vertex_update) ) {

        // Perform delta gain updates
        updateNeighbors(hn, from, to);

        // Remove all vertices that became internal from the PQ
        _border_vertices.doForAllVerticesThatBecameInternalVertices(
          [&](const HypernodeID hn) {
            ASSERT(!_border_vertices.isBorderNode(hn));
            ASSERT(_vertex_state[hn] == VertexState::ACTIVE);
            ASSERT(_pq.contains(hn));
            _pq.remove(hn, 1 - _phg.partID(hn));
            _vertex_state[hn] = VertexState::INACTIVE;
          }
        );

        // Insert all new border vertices into PQ
        _border_vertices.doForAllVerticesThatBecameBorderVertices(
          [&](const HypernodeID hn) {
          ASSERT(_border_vertices.isBorderNode(hn));
          ASSERT(_vertex_state[hn] == VertexState::INACTIVE);
          activate(hn);
        });

        performed_moves.push_back(hn);
        DBG << "Moved hypernode" << hn << "from block" << from << "to block" << to << "with gain" << gain;
        current_cut -= gain;
        current_imbalance = metrics::imbalance(_phg, _context);
        stopping_rule.update(gain);

        const bool improved_cut_within_balance = (current_cut < best_metrics.cut) &&
                                                 ( _phg.partWeight(0)
                                                   <= _context.partition.max_part_weights[0]) &&
                                                 ( _phg.partWeight(1)
                                                   <= _context.partition.max_part_weights[1]);
        const bool improved_balance_less_equal_cut = (current_imbalance < best_metrics.imbalance) &&
                                                    (current_cut <= best_metrics.cut);
        const bool move_is_feasible = ( _phg.partWeight(from) > 0) &&
                                      ( improved_cut_within_balance ||
                                        improved_balance_less_equal_cut );
        if ( move_is_feasible ) {
          DBG << GREEN << "2Way FM improved cut from" << best_metrics.cut << "to" << current_cut
              << "(Imbalance:" << current_imbalance << ")" << END;
          stopping_rule.reset();
          best_metrics.cut = current_cut;
          best_metrics.km1 = current_cut;
          best_metrics.imbalance = current_imbalance;
          min_cut_idx = performed_moves.size();
        } else {
          DBG << RED << "2Way FM decreased cut to" << current_cut
              << "(Imbalance:" << current_imbalance << ")" << END;
        }
      }
    }

    // Perform rollback to best partition found during local search
    rollback(performed_moves, min_cut_idx);

    HEAVY_REFINEMENT_ASSERT(best_metrics.cut == metrics::hyperedgeCut(_phg, false));
    HEAVY_REFINEMENT_ASSERT(best_metrics.imbalance == metrics::imbalance(_phg, _context),
           V(best_metrics.imbalance) << V(metrics::imbalance(_phg, _context)));
    return min_cut_idx > 0;
  }

 private:

  /**
   * Performs delta gain update on all non locked hyperedges and
   * state transition of hyperedges.
   */
  void updateNeighbors(const HypernodeID hn,
                       const PartitionID from,
                       const PartitionID to) {
    ASSERT(_phg.partID(hn) == to);

    for ( const HyperedgeID& he : _phg.incidentEdges(hn) ) {
      const PartitionID he_state = _he_state[he];
      if ( _phg.edgeSize(he) > 1 && he_state != HEState::LOCKED ) {
        deltaGainUpdate(he, from, to);
        // State Transition of hyperedge
        if ( he_state == HEState::FREE ) {
          // Vertex hn is the first vertex changed its block
          // in hyperedge he => free -> loose
          _he_state[he] = to;
        } else if ( he_state == from ) {
          // An other vertex already changed its block in opposite direction
          // => hyperedge he can not be removed from cut any more and therefore
          // it can not affect the gains of its pins => loose -> locked
          _he_state[he] = HEState::LOCKED;
        }
      }
    }
  }

  // ! Delta-Gain Update as decribed in [ParMar06].
  void deltaGainUpdate(const HyperedgeID he,
                       const PartitionID from,
                       const PartitionID to) {
    const HypernodeID pin_count_from_part_after_move = _phg.pinCountInPart(he, from);
    const HypernodeID pin_count_to_part_after_move = _phg.pinCountInPart(he, to);

    const bool he_became_cut_he = pin_count_to_part_after_move == 1;
    const bool he_became_internal_he = pin_count_from_part_after_move == 0;
    const bool increase_necessary = pin_count_from_part_after_move == 1;
    const bool decrease_necessary = pin_count_to_part_after_move == 2;

    if ( he_became_cut_he || he_became_internal_he ||
         increase_necessary || decrease_necessary ) {
      ASSERT(_phg.edgeSize(he) != 1, V(he));
      const HyperedgeWeight he_weight = _phg.edgeWeight(he);

      if (_phg.edgeSize(he) == 2) {
        for (const HypernodeID& pin : _phg.pins(he)) {
          if ( _vertex_state[pin] == VertexState::ACTIVE ) {
            const char factor = (_phg.partID(pin) == from ? 2 : -2);
            updatePin(pin, factor * he_weight);
          }
        }
      } else if (he_became_cut_he) {
        for (const HypernodeID& pin : _phg.pins(he)) {
          if ( _vertex_state[pin] == VertexState::ACTIVE ) {
            updatePin(pin, he_weight);
          }
        }
      } else if (he_became_internal_he) {
        for (const HypernodeID& pin : _phg.pins(he)) {
          if ( _vertex_state[pin] == VertexState::ACTIVE ) {
            updatePin(pin, -he_weight);
          }
        }
      } else {
        if ( increase_necessary || decrease_necessary ) {
          for (const HypernodeID& pin : _phg.pins(he)) {
            if ( _phg.partID(pin) == from ) {
              if ( increase_necessary && _vertex_state[pin] == VertexState::ACTIVE ) {
                updatePin(pin, he_weight);
              }
            } else if ( decrease_necessary && _vertex_state[pin] == VertexState::ACTIVE  ) {
              updatePin(pin, -he_weight);
            }
          }
        }
      }
    }
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void updatePin(const HypernodeID pin, const Gain delta) {
    const PartitionID to = 1 - _phg.partID(pin);
    ASSERT(_vertex_state[pin] == VertexState::ACTIVE, V(pin));
    ASSERT(_pq.contains(pin, to), V(pin) << V(to));
    _pq.updateKeyBy(pin, to, delta);
  }

  void activate(const HypernodeID hn) {
    if ( _border_vertices.isBorderNode(hn) ) {
      ASSERT(_vertex_state[hn] == VertexState::INACTIVE);
      const PartitionID from = _phg.partID(hn);
      const PartitionID to = 1 - from;

      ASSERT(!_pq.contains(hn, to), V(hn));
      _vertex_state[hn] = VertexState::ACTIVE;
      _pq.insert(hn, to, computeGain(hn, from, to));
      if ( _phg.partWeight(to) < _context.partition.max_part_weights[to] ) {
        _pq.enablePart(to);
      }
    }
  }

  Gain computeGain(const HypernodeID hn, const PartitionID from, const PartitionID to) {
    ASSERT(_phg.partID(hn) == from);
    ASSERT(1 - from == to);
    Gain gain = 0;
    for ( const HyperedgeID& he : _phg.incidentEdges(hn) ) {
      if ( _phg.edgeSize(he) > 1 ) {
        if ( _phg.pinCountInPart(he, to) == 0 ) {
          gain -= _phg.edgeWeight(he);
        }
        if ( _phg.pinCountInPart(he, from) == 1 ) {
          gain += _phg.edgeWeight(he);
        }
      }
    }
    return gain;
  }

  void rollback(const parallel::scalable_vector<HypernodeID>& performed_moves,
                const size_t min_cut_idx) {
    for ( size_t i = min_cut_idx; i < performed_moves.size(); ++i ) {
      const HypernodeID hn = performed_moves[i];
      const PartitionID from = _phg.partID(hn);
      const PartitionID to = 1 - from;
      _phg.changeNodePart(hn, from, to);
    }
  }

  bool verifyPQState() const {
    for ( const HypernodeID& hn : _phg.nodes() ) {
      const PartitionID to = 1 - _phg.partID(hn);
      if ( _border_vertices.isBorderNode(hn) && _vertex_state[hn] != VertexState::MOVED ) {
        if ( !_pq.contains(hn, to) ) {
          LOG << "Hypernode" << hn << "is a border node and should be contained in the PQ";
          return false;
        }
        if ( _vertex_state[hn] != VertexState::ACTIVE ) {
          LOG << "Hypernode" << hn << "is a border node and its not moved and its state should be ACTIVE";
          return false;
        }
      } else if ( !_border_vertices.isBorderNode(hn) && _vertex_state[hn] != VertexState::MOVED ) {
        if ( _pq.contains(hn, to) ) {
          LOG << "Hypernode" << hn << "is not a border node and should be not contained in PQ";
          return false;
        }
        if ( _vertex_state[hn] != VertexState::INACTIVE ) {
          LOG << "Hypernode" << hn << "is not a border node and its not moved and its state should be INACTIVE";
          return false;
        }
      }
    }
    return true;
  }

  PartitionedHypergraph& _phg;
  const Context& _context;

  KWayPriorityQueue _pq;
  BorderVertexTracker _border_vertices;
  parallel::scalable_vector<VertexState> _vertex_state;
  parallel::scalable_vector<PartitionID> _he_state;
};

}

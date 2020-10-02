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

#include "mt-kahypar/partition/refinement/fm/sequential_twoway_fm_refiner.h"

#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/fm/stop_rule.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {

bool SequentialTwoWayFmRefiner::refine(kahypar::Metrics& best_metrics) {

  // Activate all border nodes
  _pq.clear();
  _border_vertices.initialize(_phg);
  utils::Randomize::instance().shuffleVector(_nodes, sched_getcpu());
  for ( const HypernodeID& hn : _nodes ) {
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
          _context.partition.max_part_weights[to], []{}, border_vertex_update) ) {

      // Perform delta gain updates
      updateNeighbors(hn, from, to);
      updatePQState(from, to);

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

void SequentialTwoWayFmRefiner::activate(const HypernodeID hn) {
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

/**
 * Performs delta gain update on all non locked hyperedges and
 * state transition of hyperedges.
 */
void SequentialTwoWayFmRefiner::updateNeighbors(const HypernodeID hn,
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
void SequentialTwoWayFmRefiner::deltaGainUpdate(const HyperedgeID he,
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

void SequentialTwoWayFmRefiner::updatePin(const HypernodeID pin, const Gain delta) {
  const PartitionID to = 1 - _phg.partID(pin);
  ASSERT(_vertex_state[pin] == VertexState::ACTIVE, V(pin));
  ASSERT(_pq.contains(pin, to), V(pin) << V(to));
  _pq.updateKeyBy(pin, to, delta);
}

void SequentialTwoWayFmRefiner::updatePQState(const PartitionID from,
                                              const PartitionID to) {
  if (_phg.partWeight(to) >= _context.partition.max_part_weights[to] ) {
    _pq.disablePart(to);
  }
  if (_phg.partWeight(from) < _context.partition.max_part_weights[from] ) {
    _pq.enablePart(from);
  }
}

Gain SequentialTwoWayFmRefiner::computeGain(const HypernodeID hn, const PartitionID from, const PartitionID to) {
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

void SequentialTwoWayFmRefiner::rollback(const parallel::scalable_vector<HypernodeID>& performed_moves,
              const size_t min_cut_idx) {
  for ( size_t i = min_cut_idx; i < performed_moves.size(); ++i ) {
    const HypernodeID hn = performed_moves[i];
    const PartitionID from = _phg.partID(hn);
    const PartitionID to = 1 - from;
    _phg.changeNodePart(hn, from, to);
  }
}

bool SequentialTwoWayFmRefiner::verifyPQState() const {
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

} // namespace mt_kahypar

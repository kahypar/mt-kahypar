/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2018 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2018 Tobias Heuer <tobias.heuer@live.com>
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

#include <vector>

#include "tbb/enumerable_thread_specific.h"

#include "kahypar/meta/policy_registry.h"
#include "kahypar/meta/typelist.h"
#include "kahypar/partition/metrics.h"
#include "kahypar/partition/context_enum_classes.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {


template <class Derived = Mandatory,
          class HyperGraph = Mandatory>
class GainPolicy : public kahypar::meta::PolicyBase {

using DeltaGain = tbb::enumerable_thread_specific<Gain>;
using TmpScores = tbb::enumerable_thread_specific<parallel::scalable_vector<Gain>>;

 public:
  GainPolicy(HyperGraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _context(context),
    _deltas(0),
    _tmp_scores(context.partition.k, 0) { }

  Move computeMaxGainMove(const HypernodeID hn) {
    return static_cast<Derived*>(this)->computeMaxGainMoveImpl(hn);
  }

  inline void computeDeltaForHyperedge(const HyperedgeWeight edge_weight,
                                       const HypernodeID edge_size,
                                       const HypernodeID pin_count_in_from_part_after,
                                       const HypernodeID pin_count_in_to_part_after) {
    static_cast<Derived*>(this)->computeDeltaForHyperedgeImpl(
      edge_weight, edge_size, pin_count_in_from_part_after, pin_count_in_to_part_after);
  }

  // ! Returns the delta in the objective function for all moves
  // ! performed by the calling thread relative to the last call
  // ! reset()
  Gain localDelta() {
    return _deltas.local();
  }

  // ! Returns the overall delta of all moves performed by
  // ! all threads relative to the last call of reset()
  Gain delta() const {
    Gain overall_delta = 0;
    for ( const Gain& delta : _deltas ) {
      overall_delta += delta;
    }
    return overall_delta;
  }

  void reset() {
    for ( Gain& delta : _deltas ) {
      delta = 0;
    }
  }

 protected:
  HyperGraph& _hg;
  const Context& _context;
  DeltaGain _deltas;
  TmpScores _tmp_scores;
};

template <class HyperGraph = Mandatory>
class Km1Policy : public GainPolicy<Km1Policy<HyperGraph>, HyperGraph> {

using Base = GainPolicy<Km1Policy<HyperGraph>, HyperGraph>;

static constexpr bool enable_heavy_assert = false;

 public:
  Km1Policy(HyperGraph& hypergraph,
            const Context& context,
            bool disable_randomization = false) :
    Base(hypergraph, context),
    _disable_randomization(disable_randomization) { }

  Move computeMaxGainMoveImpl(const HypernodeID hn) {
    HEAVY_REFINEMENT_ASSERT([&] {
      for (PartitionID k = 0; k < _context.partition.k; ++k) {
        if ( _tmp_scores.local()[k] != 0 ) {
          return false;
        }
      }
      return true;
    }(), "Scores and valid parts not correctly reset");

    PartitionID from = _hg.partID(hn);
    parallel::scalable_vector<Gain>& tmp_scores = _tmp_scores.local();
    for ( const HyperedgeID& he : _hg.incidentEdges(hn) ) {
      HypernodeID pin_count_in_from_part = _hg.pinCountInPart(he, from);
      HyperedgeWeight he_weight = _hg.edgeWeight(he);
      if ( pin_count_in_from_part == 1 ) {
        // In case, there is only one pin left in the from part,
        // we can decrease the connectivity metric by moving the
        // only pin left (vertex 'hn') to one of its incident blocks
        // also contained in the hyperedge
        for ( const PartitionID& to : _hg.connectivitySet(he) ) {
          if ( from != to ) {
            tmp_scores[to] -= he_weight;
          }
        }
      } else {
        // In case, there are more than one pin left in from part, than
        // we would increase the connectivity, if we would move vertex hn
        // to one block that is not contained in the hyperedge.
        for ( PartitionID to = 0; to < _context.partition.k; ++to ) {
          if ( from != to && _hg.pinCountInPart(he, to) == 0 ) {
            tmp_scores[to] += he_weight;
          }
        }
      }
    }

    Move best_move { from, from, 0 };
    HypernodeWeight hn_weight = _hg.nodeWeight(hn);
    int cpu_id = sched_getcpu();
    utils::Randomize& rand = utils::Randomize::instance();
    for ( PartitionID to = 0; to < _context.partition.k; ++to ) {
      if ( from != to ) {
        bool new_best_gain = ( tmp_scores[to] < best_move.gain ) ||
                            ( tmp_scores[to] == best_move.gain &&
                              !_disable_randomization &&
                              rand.flipCoin(cpu_id) );
        if ( new_best_gain && _hg.localPartWeight(to) + hn_weight <=
             _context.partition.max_part_weights[to] ) {
          best_move.to = to;
          best_move.gain = tmp_scores[to];
        }
      }
      tmp_scores[to] = 0;
    }
    return best_move;
  }

  inline void computeDeltaForHyperedgeImpl(const HyperedgeWeight edge_weight,
                                           const HypernodeID edge_size,
                                           const HypernodeID pin_count_in_from_part_after,
                                           const HypernodeID pin_count_in_to_part_after) {
    _deltas.local() += HyperGraph::km1Delta(edge_weight, edge_size,
      pin_count_in_from_part_after, pin_count_in_to_part_after);
  }

  using Base::_hg;
  using Base::_context;
  using Base::_deltas;
  using Base::_tmp_scores;
  bool _disable_randomization;
};


template <class HyperGraph = Mandatory>
class CutPolicy : public GainPolicy<CutPolicy<HyperGraph>, HyperGraph> {

using Base = GainPolicy<CutPolicy<HyperGraph>, HyperGraph>;

static constexpr bool enable_heavy_assert = false;

 public:
  CutPolicy(HyperGraph& hypergraph,
            const Context& context,
            bool disable_randomization = false) :
    Base(hypergraph, context),
    _disable_randomization(disable_randomization) { }

  Move computeMaxGainMoveImpl(const HypernodeID hn) {
    HEAVY_REFINEMENT_ASSERT([&] {
      for (PartitionID k = 0; k < _context.partition.k; ++k) {
        if ( _tmp_scores.local()[k] != 0 ) {
          return false;
        }
      }
      return true;
    }(), "Scores and valid parts not correctly reset");

    PartitionID from = _hg.partID(hn);
    parallel::scalable_vector<Gain>& tmp_scores = _tmp_scores.local();
    Gain internal_weight = 0;
    for ( const HyperedgeID& he : _hg.incidentEdges(hn) ) {
      PartitionID connectivity = _hg.connectivity(he);
      HypernodeID pin_count_in_from_part = _hg.pinCountInPart(he, from);
      HyperedgeWeight weight = _hg.edgeWeight(he);
      if ( connectivity == 1 ) {
        ASSERT(_hg.edgeSize(he) > 1);
        // In case, the hyperedge is a non-cut hyperedge, we would increase
        // the cut, if we move vertex hn to an other block.
        internal_weight += weight;
      } else if ( connectivity == 2 && pin_count_in_from_part == 1 ) {
        for ( const PartitionID& to : _hg.connectivitySet(he) ) {
          // In case there are only two blocks contained in the current
          // hyperedge and only one pin left in the from part of the hyperedge,
          // we would make the current hyperedge a non-cut hyperedge when moving
          // vertex hn to the other block.
          if ( from != to ) {
            tmp_scores[to] -= weight;
          }
        }
      }
    }

    Move best_move { from, from, 0 };
    HypernodeWeight hn_weight = _hg.nodeWeight(hn);
    int cpu_id = sched_getcpu();
    utils::Randomize& rand = utils::Randomize::instance();
    for ( PartitionID to = 0; to < _context.partition.k; ++to ) {
      if ( from != to ) {
        Gain score = tmp_scores[to] + internal_weight;
        bool new_best_gain = ( score < best_move.gain ) ||
                             ( score == best_move.gain &&
                               !_disable_randomization &&
                               rand.flipCoin(cpu_id) );
        if ( new_best_gain && _hg.localPartWeight(to) + hn_weight <=
             _context.partition.max_part_weights[to] ) {
          best_move.to = to;
          best_move.gain = score;
        }
      }
      tmp_scores[to] = 0;
    }
    return best_move;
  }

  inline void computeDeltaForHyperedgeImpl(const HyperedgeWeight edge_weight,
                                           const HypernodeID edge_size,
                                           const HypernodeID pin_count_in_from_part_after,
                                           const HypernodeID pin_count_in_to_part_after) {
    _deltas.local() += HyperGraph::cutDelta(edge_weight, edge_size,
      pin_count_in_from_part_after, pin_count_in_to_part_after);
  }

  using Base::_hg;
  using Base::_context;
  using Base::_deltas;
  using Base::_tmp_scores;
  bool _disable_randomization;
};

}  // namespace mt_kahypar

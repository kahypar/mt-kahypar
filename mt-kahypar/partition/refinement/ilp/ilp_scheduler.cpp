/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019, 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "ilp_scheduler.h"

#include "tbb/parallel_sort.h"

#include "mt-kahypar/partition/refinement/policies/gain_policy.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/parallel/stl/scalable_queue.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/progress_bar.h"

namespace mt_kahypar {

bool ILPScheduler::refine() {

  utils::Timer::instance().start_timer("compute_gains", "Compute Gains");
  const bool compute_gains = ( _context.refinement.ilp.vertex_selection_strategy ==
                               ILPVertexSelectionStrategy::gain ||
                               _context.refinement.ilp.vertex_selection_strategy ==
                               ILPVertexSelectionStrategy::top_vertices );
  _gains.assign(_phg.initialNumNodes(), VertexGain { kInvalidHypernode, kInvalidGain, false });
  if ( compute_gains ) {
    // Compute gains of each vertex
    computeGains();
  } else {
    // Only determine which vertices are border vertices
    _phg.doParallelForAllNodes([&](const HypernodeID& hn) {
      _gains[hn].hn = hn;
      _gains[hn].is_border_vertex = _phg.isBorderNode(hn);
    });
  }

  if ( _context.refinement.ilp.vertex_selection_strategy != ILPVertexSelectionStrategy::localized ) {
    tbb::parallel_sort(_gains.begin(), _gains.end(),
      [&](const VertexGain& lhs, const VertexGain& rhs) {
        return ( lhs.is_border_vertex > rhs.is_border_vertex ||
              ( lhs.is_border_vertex == rhs.is_border_vertex && lhs.gain > rhs.gain ) ||
              ( lhs.is_border_vertex == rhs.is_border_vertex && lhs.gain == rhs.gain && lhs.hn < rhs.hn ) );
      });
  } else {
    utils::Randomize::instance().parallelShuffleVector(_gains, 0UL, _gains.size());
  }
  utils::Timer::instance().stop_timer("compute_gains");

  _visited_hns.assign(_phg.initialNumNodes(), false);
  _visited_hes = kahypar::ds::FastResetFlagArray<>(_phg.initialNumEdges());
  if ( _context.refinement.ilp.vertex_selection_strategy != ILPVertexSelectionStrategy::localized ) {
    utils::Timer::instance().start_timer("select_vertices", "Select ILP Vertices");
    _num_vertices = 0;
    _num_hyperedges = 0;
    _num_pins = 0;
    _k = 0;
    vec<HypernodeID> nodes;
    if ( _context.refinement.ilp.vertex_selection_strategy == ILPVertexSelectionStrategy::boundary ) {
      bfs(nodes, 0, _gains.size(), kInvalidGain, std::numeric_limits<int>::max());
    } else if ( _context.refinement.ilp.vertex_selection_strategy == ILPVertexSelectionStrategy::gain ) {
      size_t pos = bfs(nodes, 0, _gains.size(), _context.refinement.ilp.min_gain, std::numeric_limits<int>::max());
      if ( estimatedNumberOfNonZeros() < _context.refinement.ilp.max_non_zeros ) {
        bfs(nodes, pos, _gains.size(), kInvalidGain, std::numeric_limits<int>::max());
      }
    } else if ( _context.refinement.ilp.vertex_selection_strategy == ILPVertexSelectionStrategy::top_vertices ) {
      size_t idx = 0;
      while ( estimatedNumberOfNonZeros() < _context.refinement.ilp.max_non_zeros &&
              idx < _gains.size() && _gains[idx].is_border_vertex ) {
        bfs(nodes, idx, idx + 1, kInvalidGain, _context.refinement.ilp.max_bfs_distance);
        ++idx;
      }
    }
    utils::Timer::instance().stop_timer("select_vertices");

    // Solve ILP
    _solver.solve(nodes);
  } else {
    HyperedgeWeight initial_objective = metrics::objective(_phg, _context.partition.objective);
    HyperedgeWeight current_objective = initial_objective;
    HyperedgeWeight max_allowed_objective = initial_objective;
    double current_epsilon = 0.0;
    while ( maxPartWeight() > _context.partition.perfect_balance_part_weights[0] ) {
      _phg.doParallelForAllNodes([&](const HypernodeID& hn) {
        _gains[hn].hn = hn;
        _gains[hn].is_border_vertex = _phg.isBorderNode(hn);
      });
      utils::Randomize::instance().parallelShuffleVector(_gains, 0UL, _gains.size());

      _marked_hns.assign(_phg.initialNumNodes(), false);
      _visited_hns.assign(_phg.initialNumNodes(), false);
      utils::ProgressBar ilp_progress(_phg.initialNumNodes(),
                                      _context.refinement.ilp.minimize_balance ?
                                      maxPartWeight() : current_objective,
                                      _context.partition.verbose_output &&
                                      _context.partition.enable_progress_bar && !debug);
      for ( size_t i = 0; i < _gains.size(); ++i ) {
        const HypernodeID hn = _gains[i].hn;
        if ( _gains[i].is_border_vertex && !_marked_hns[hn] ) {
          utils::Timer::instance().start_timer("select_vertices", "Select ILP Vertices");
          _num_vertices = 0;
          _num_hyperedges = 0;
          _num_pins = 0;
          _k = 0;
          _visited_hes.reset();
          vec<HypernodeID> nodes;
          bfs(nodes, hn, _context.refinement.ilp.max_bfs_distance);
          utils::Timer::instance().stop_timer("select_vertices");

          // Solve ILP
          if ( estimatedNumberOfNonZeros() >= _context.refinement.ilp.min_non_zeros ) {
            const HypernodeWeight max_part_weight_before = maxPartWeight();
            HyperedgeWeight delta = _solver.solve(
              nodes, !debug /* supress output */, max_allowed_objective - current_objective);
            const HypernodeWeight max_part_weight_after = maxPartWeight();
            current_objective -= delta;
            ilp_progress.setObjective(
              _context.refinement.ilp.minimize_balance ? max_part_weight_after : current_objective);
            ilp_progress += nodes.size();
            if ( _context.refinement.ilp.minimize_balance &&
                max_part_weight_before == max_part_weight_after ) {
              current_epsilon += 0.01;
              max_allowed_objective = (1.0 + current_epsilon) * initial_objective;
            }
            // LOG << V(initial_objective) << V(current_objective) << V(max_allowed_objective) << V(max_part_weight_after);
          }
        }

        if ( maxPartWeight() == _context.partition.perfect_balance_part_weights[0] ) {
          break;
        }
      }
      ilp_progress += (_phg.initialNumNodes() - ilp_progress.count());

      if ( !_context.refinement.ilp.minimize_balance ) {
        break;
      }
    }
  }

  return true;
}

void ILPScheduler::computeGainsForConnectivityMetric() {
  Km1Policy<PartitionedHypergraph> _gain_calculator(_context);
  _phg.doParallelForAllNodes([&](const HypernodeID& hn) {
    _gains[hn].hn = hn;
    _gains[hn].is_border_vertex = _phg.isBorderNode(hn);
    if ( _gains[hn].is_border_vertex ) {
      Move move = _gain_calculator.computeMaxGainMove(_phg, hn, false, true);
      _gains[hn].gain = -move.gain;
    }
  });
}

void ILPScheduler::computeGainsForCutMetric() {
  CutPolicy<PartitionedHypergraph> _gain_calculator(_context);
  _phg.doParallelForAllNodes([&](const HypernodeID& hn) {
    _gains[hn].hn = hn;
    _gains[hn].is_border_vertex = _phg.isBorderNode(hn);
    if ( _gains[hn].is_border_vertex ) {
      Move move = _gain_calculator.computeMaxGainMove(_phg, hn, false, true);
      _gains[hn].gain = -move.gain;
    }
  });
}

// ! Returns the number of vertices initially inserted into the bfs queue
size_t ILPScheduler::bfs(vec<HypernodeID>& nodes,
                         const size_t gains_start_idx,
                         const size_t gains_end_idx,
                         const Gain min_gain,
                         const int max_distance) {
  ASSERT(gains_start_idx < gains_end_idx);
  ASSERT(gains_end_idx <= _gains.size());
  vec<bool> visited_k(_context.partition.k, false);
  parallel::scalable_queue<HypernodeID> q;
  parallel::scalable_queue<HypernodeID> next_q;

  // Insert vertices into BFS queue
  size_t idx = gains_start_idx;
  for ( ; idx < gains_end_idx; ++idx ) {
    const HypernodeID hn = _gains[idx].hn;
    if ( _gains[idx].is_border_vertex && _gains[idx].gain >= min_gain && !_visited_hns[hn] ) {
      q.push(hn);
      _visited_hns[hn] = true;
    } else {
      break;
    }
  }

  // Perform BFS
  int current_distance = 0;
  while ( !q.empty() ) {
    const HypernodeID hn = q.front();
    q.pop();

    nodes.push_back(hn);
    ++_num_vertices;
    _num_pins += _phg.nodeDegree(hn);

    for ( const HyperedgeID& he : _phg.incidentEdges(hn) ) {
      if ( !_visited_hes[he] ) {
        for ( const HypernodeID& pin : _phg.pins(he) ) {
          if ( !_visited_hns[pin] ) {
            next_q.push(pin);
            _visited_hns[pin] = true;
          }
        }
        for ( const PartitionID& i : _phg.connectivitySet(he) ) {
          if ( !visited_k[i] ) {
            ++_k;
            visited_k[i] = true;
          }
        }
        ++_num_hyperedges;
        _visited_hes.set(he, true);
      }
    }

    if ( estimatedNumberOfNonZeros() >= _context.refinement.ilp.max_non_zeros ) {
      break;
    }

    if ( q.empty() ) {
      ++current_distance;
      if ( current_distance <= max_distance ) {
        q.swap(next_q);
      }
    }
  }

  return idx;
}

namespace {
  struct QueueElement {
    HypernodeID hn;
    int distance;
  };
} // namespace

// ! Returns the number of vertices initially inserted into the bfs queue
void ILPScheduler::bfs(vec<HypernodeID>& nodes,
                       const HypernodeID start_hn,
                       const int max_distance) {
  vec<bool> visited_k(_context.partition.k, false);
  parallel::scalable_queue<QueueElement> q;
  q.push(QueueElement { start_hn, 0 });
  _visited_hns[start_hn] = true;

  while ( !q.empty() ) {
    const HypernodeID hn = q.front().hn;
    const int current_distance = _phg.isBorderNode(hn) ? 0 : q.front().distance;
    q.pop();

    nodes.push_back(hn);
    _marked_hns[hn] = true;
    ++_num_vertices;
    _num_pins += _phg.nodeDegree(hn);

    if ( current_distance < max_distance ) {
      for ( const HyperedgeID& he : _phg.incidentEdges(hn) ) {
        if ( !_visited_hes[he] ) {
          for ( const HypernodeID& pin : _phg.pins(he) ) {
            if ( !_visited_hns[pin] ) {
              q.push(QueueElement { pin, current_distance + 1 });
              _visited_hns[pin] = true;
            }
          }
          for ( const PartitionID& i : _phg.connectivitySet(he) ) {
            if ( !visited_k[i] ) {
              ++_k;
              visited_k[i] = true;
            }
          }
          ++_num_hyperedges;
          _visited_hes.set(he, true);
        }
      }
    }

    if ( estimatedNumberOfNonZeros() >= _context.refinement.ilp.max_non_zeros ) {
      break;
    }
  }

  // Enable all vertices not added to ILP again
  while ( !q.empty() ) {
    const HypernodeID hn = q.front().hn;
    q.pop();
    if ( !_marked_hns[hn] ) {
      _visited_hns[hn] = false;
    }
  }
}


} // namespace mt_kahypar
/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
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

#include "mt-kahypar/partition/refinement/rebalancing/md_rebalancer.h"

#include <tbb/parallel_sort.h>
#include <ranges>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/cast.h"

namespace mt_kahypar{

  template <typename GraphAndGainTypes>
  bool MDRebalancer<GraphAndGainTypes>::refineInternal(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                       vec<vec<Move>>* moves_by_part,
                                                       vec<Move>* moves_linear,
                                                       Metrics& best_metrics) {
                                                      
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
    resizeDataStructuresForCurrentK();
    std::cout << "teststart\n" << phg.initialNumNodes() << "\n";
    // This function is passed as lambda to the changeNodePart function and used
    // to calculate the "real" delta of a move (in terms of the used objective function).
    auto objective_delta = [&](const SynchronizedEdgeUpdate& sync_update) {
      _gain.computeDeltaForHyperedge(sync_update);
    };

    // If partition is imbalanced, rebalancer is activated
    if ( !metrics::isBalanced(phg, _context) ) {
      DBG << "Starting multi-dimensional rebalancer";  // only printed if debug=true in header
      _gain.reset();

    auto balance_gain = [&](const PartitionedHypergraph& phg, HypernodeID node, PartitionID from, PartitionID to){
      double gain = 0;
      for(int i = 0; i < dimension; i++){
        gain += std::max(0, std::min(phg.nodeWeight(node).weights[i], phg.partWeight(to).weights[i] + 
        phg.nodeWeight(node).weights[i] - _context.partition.max_part_weights[to].weights[i])) * _context.partition.max_part_weights_inv[to][i]
        - std::max(0, std::min(phg.nodeWeight(node).weights[i], phg.partWeight(from).weights[i] - _context.partition.max_part_weights[from].weights[i])) * _context.partition.max_part_weights_inv[from][i];
      }
      return gain;
    };
    std::cout << "start nodes\n";
    std::vector<tbb::concurrent_vector<std::pair<HypernodeID, HypernodeID>>> nodes_sorted;
    std::vector<HypernodeWeight> exceed;
    exceed.resize(phg.k());
    auto calcExceed = [&](PartitionID p){
      for(int j = 0; j < dimension; j++){
        exceed[p].weights[j] = phg.partWeight(p).weights[j] - _context.partition.max_part_weights[p].weights[j];
      }
    };
    for(PartitionID p = 0; p < phg.k(); p++){
      calcExceed(p);
    }
    std::vector<std::vector<HypernodeID>> indices;
    /*std::vector<std::vector<std::list<HypernodeID>>> partition_lists;*/
    nodes_sorted.resize(mt_kahypar::dimension);
    indices.resize(mt_kahypar::dimension);
    /*partition_lists.resize(mt_kahypar::dimension);*/
    for(int i = 0; i < dimension; i++){
      indices[i].resize(phg.k());
      /*partition_lists[i].resize(phg.k());*/
    }
    phg.doParallelForAllNodes([&](const HypernodeID hn){
      for(int i = 0; i < dimension; i++){
        ASSERT(hn < phg.initialNumNodes());
        nodes_sorted[i].push_back({phg.nodeWeight(hn).weights[i], hn});
      }
    });
    for(int i = 0; i < dimension; i++){
      ASSERT(nodes_sorted[i].size() == phg.initialNumNodes());
    }
    tbb::parallel_for(UL(0), mt_kahypar::dimension, [&](const size_t i){tbb::parallel_sort(nodes_sorted[i].begin(), nodes_sorted[i].end());});
    auto getBorder = [&](int dimension, HypernodeWeight nw, HypernodeID offset = 0){
      HypernodeID left = offset;
      HypernodeID right = nodes_sorted[dimension].size();
      while(left < right - 1){
        HypernodeID middle = (left + right) / 2;
      if(nodes_sorted[dimension][middle].first <= nw.weights[dimension]){
        left = middle;
      }
      else{
        right = middle;
      }
      }
      return right;
      /*return std::upper_bound(nodes_sorted[dimension].begin() + offset, nodes_sorted[dimension].end(), 
        {std::abs(nw.weights[dimension]), std::numeric_limits<HypernodeID>::min()})
         - nodes_sorted[dimension].begin();*/
    };
    tbb::parallel_for(UL(0), mt_kahypar::dimension, [&](const size_t i){
      tbb::parallel_for(UL(0), (size_t)phg.k(), [&](const size_t j){
        indices[i][j] = getBorder(i, exceed[j]);
      });     
    });

    std::cout << "start queue\n";
    MoveQueue queue;
    queue.initialize(phg.initialNumNodes(), phg.k());
    phg.doParallelForAllNodes([&](const HypernodeID hn){
      std::vector<Move_md> moves = _gain.allGains(phg, hn);
      for(int i = 0; i < moves.size(); i++){
        ASSERT(phg.partID(hn) != moves[i].to);
        ASSERT(phg.partID(hn) != -1);
        queue.insert_without_updating({hn, {moves[i].to, {moves[i].gain_and_balance, moves[i].gain, moves[i].balance}}});
      }
    });
    queue.check();
    std::cout << "end queue\n";
    PartitionID imbalanced = 0;
    for(PartitionID p = 0; p < phg.k(); p++){
      if(phg.partWeight(p) > _context.partition.max_part_weights[p]){
        imbalanced++;
      }
    }
    auto add = [&](HypernodeID hn, PartitionID p, int dimension, int32_t balance){
      queue.addBalance(hn, p, static_cast<double>(balance) * _context.partition.max_part_weights_inv[phg.partID(hn)][dimension]);
    };
    auto addAll = [&](HypernodeID hn, int dimension, int32_t balance){
      for(PartitionID p = 0; p < phg.k(); p++){
        if(p != phg.partID(hn)){
          queue.addBalance(hn, p, static_cast<double>(balance) * _context.partition.max_part_weights_inv[phg.partID(hn)][dimension]);
        }
      }
    };

    while(imbalanced != 0 && !queue.isEmpty()){
      queue.top_moves.swap_down_counter = 0;
      queue.top_moves.swap_up_counter = 0;
      std::pair<HypernodeID, std::pair<PartitionID, Move_internal>> max_move = queue.deleteMax();
      ASSERT(max_move.first < phg.initialNumNodes());
      ASSERT(phg.partID(max_move.first) != -1);
      ASSERT(max_move.second.first != -1);
      HypernodeID node = max_move.first;
      PartitionID from = phg.partID(node);
      PartitionID to = max_move.second.first;
      double balance = balance_gain(phg, node, from, to);
      if(balance > max_move.second.second.balance + 0.000000000000000000001){
        queue.changeBalance({node, {to, balance}});
        queue.update(node);
        continue;
      }
      Move_internal gains = max_move.second.second;
      HypernodeWeight from_old = phg.partWeight(from);
      HypernodeWeight to_old = phg.partWeight(to);
      HypernodeWeight exceed_from_old = exceed[from];
      HypernodeWeight exceed_to_old = exceed[to];

      Move move = {phg.partID(node), to, node, 0};
      imbalanced += (phg.partWeight(phg.partID(node)) - phg.nodeWeight(node) > _context.partition.max_part_weights[phg.partID(node)])
        - (phg.partWeight(phg.partID(node)) > _context.partition.max_part_weights[phg.partID(node)])
        - (phg.partWeight(to) > _context.partition.max_part_weights[to])
        + (phg.partWeight(to) + phg.nodeWeight(node) > _context.partition.max_part_weights[to]);
      phg.changeNodePart(node, phg.partID(node), to, objective_delta);
      calcExceed(from);
      calcExceed(to);
      tbb::concurrent_vector<HypernodeID> changed_nodes = _gain.getChangedMoves(phg, {move.from, move.to, move.node, 0}, &queue);
      /*for(PartitionID p = 0; p < phg.k(); p++){
        queue.addToGain({node, {p, -gains.gain}});
      }*/
      for(int i = 0; i < dimension; i++){
        if(exceed[from].weights[i] > 0){
          for(HypernodeID h = nodes_sorted[i].size() - 1; h >= indices[i][from]; h--){
            HypernodeID hn = nodes_sorted[i][h].second;
            ASSERT(hn < phg.initialNumNodes());
            if(phg.partID(hn) == from){
              addAll(hn, i, phg.nodeWeight(node).weights[i]);
              changed_nodes.push_back(hn);
            }
          }
          while(indices[i][from] > 0 && nodes_sorted[i][indices[i][from] - 1].first > exceed[from].weights[i]){
            indices[i][from]--;
            HypernodeID hn = nodes_sorted[i][indices[i][from]].second;
            ASSERT(hn < phg.initialNumNodes());
            if(phg.partID(hn) == from){
              addAll(hn, i, (nodes_sorted[i][indices[i][from]].first - exceed[from].weights[i]));
              changed_nodes.push_back(hn);
            }
          }
        }
        else if(exceed_from_old.weights[i] > 0){
          std::cout << "up " << from << "\n";
          indices[i][from] = getBorder(i, exceed[from]);
          for(HypernodeID hn = 0; hn < nodes_sorted[i].size(); hn++){
            ASSERT(hn < phg.initialNumNodes());
            if(phg.partID(hn) == from){
              addAll(hn, i, std::min(phg.nodeWeight(hn).weights[i], exceed_from_old.weights[i]));
              changed_nodes.push_back(hn);
            }
            else{
              add(hn, from, i, -std::min(phg.nodeWeight(hn).weights[i], -exceed[from].weights[i]));
              changed_nodes.push_back(hn);
            }
          }
        }
        else{
          while(indices[i][from] < nodes_sorted[i].size() && nodes_sorted[i][indices[i][from]].first <= -exceed[from].weights[i]){
            HypernodeID hn = nodes_sorted[i][indices[i][from]].second;
            ASSERT(hn < phg.initialNumNodes());
            if(phg.partID(hn) != from){
              add(hn, from, i, -(nodes_sorted[i][indices[i][from]].first + exceed_from_old.weights[i]));
              changed_nodes.push_back(hn);
            }
            indices[i][from]++;
          }
          for(HypernodeID h = indices[i][from]; h < nodes_sorted[0].size(); h++){
            HypernodeID hn = nodes_sorted[i][h].second;
            ASSERT(hn < phg.initialNumNodes());
            if(phg.partID(hn) != from){
              add(hn, from, i, -nodes_sorted[i][h].first);
              changed_nodes.push_back(hn);
            }
          }
        }
        if(exceed[to].weights[i] <= 0){
          for(HypernodeID h = nodes_sorted[i].size() - 1; h >= indices[i][to];  h--){
            HypernodeID hn = nodes_sorted[i][h].second;
            ASSERT(hn < phg.initialNumNodes(), hn);
            if(phg.partID(hn) != to){
              add(hn, to, i, phg.nodeWeight(node).weights[i]);
              changed_nodes.push_back(hn);
            }
          }
          while(indices[i][to] > 0 && nodes_sorted[i][indices[i][to] - 1].first > exceed[to].weights[i]){
            indices[i][to]--;
            HypernodeID hn = nodes_sorted[i][indices[i][to]].second;
            ASSERT(hn < phg.initialNumNodes());
            if(phg.partID(hn) != to){
              add(hn, to, i, nodes_sorted[i][indices[i][to]].first + exceed[to].weights[i]);
              changed_nodes.push_back(hn);
            }
          }          
        }
        else if(exceed_to_old.weights[i] <= 0){
          std::cout << "down " << to << "\n";
          indices[i][to] = getBorder(i, exceed[to]);
          for(HypernodeID hn = 0; hn < nodes_sorted[0].size(); hn++){
            ASSERT(hn < phg.initialNumNodes());
            if(phg.partID(hn) == to){
              addAll(hn, i, -std::min(exceed[to].weights[i], phg.nodeWeight(hn).weights[i]));
              changed_nodes.push_back(hn);
            }
            else{
              add(hn, to, i, std::min(-exceed_to_old.weights[i], phg.nodeWeight(hn).weights[i]));
              changed_nodes.push_back(hn);
            }
          }
        }
        else{
          while(indices[i][to] < nodes_sorted[0].size() && nodes_sorted[i][indices[i][to]].first <= exceed[to].weights[i]){
            HypernodeID hn = nodes_sorted[i][indices[i][to]].second;
            ASSERT(hn < phg.initialNumNodes());
            if(phg.partID(hn) == to){
              addAll(hn, i, -exceed_to_old.weights[i]);
              changed_nodes.push_back(hn);
            }
            indices[i][to]++;
          }
          for(HypernodeID h = indices[i][to]; h < nodes_sorted[0].size(); h++){
            HypernodeID hn = nodes_sorted[i][h].second;
            ASSERT(hn < phg.initialNumNodes());
            if(phg.partID(hn) == to){
              addAll(hn, i, nodes_sorted[i][h].first);
              changed_nodes.push_back(hn);
            }
          }
        }

        /*if(phg.partWeight(from).weights[i] + phg.nodeWeight(node).weights[i] > _context.partition.max_part_weights[from].weights[i]){
          if(phg.partWeight(from).weights[i] > _context.partition.max_part_weights[from].weights[i]){
            for(HypernodeID j = nodes_sorted[0].size() - 1; j >= indices[i][from]; j--){
              HypernodeID hn = phg.partID(nodes_sorted[i][j].second);
              if(phg.partID(hn) == from){
                queue.checkAll(hn, i);
              }
            }
            while(indices[i][from] > 0 && nodes_sorted[i][indices[i][from] - 1] > std::abs(_context.partition.max_part_weights[from][i] - phg.partWeight(from).weights[i])){
              indices[i][from]--;
              HypernodeID hn = phg.partID(nodes_sorted[i][indices[i][from]].second);
              if(phg.partID(hn) == from){
                queue.checkAll(hn, i);
              }
            }
          }
          else{
            indices[i][from] = std::ranges::upper_bound(nodes_sorted[i].begin(), nodes_sorted[i].end(), 
              {std::abs(_context.partition.max_part_weights[from][i] - phg.partWeight(from).weights[i]) , std::numeric_limits<HypernodeID>::min()})
               - nodes_sorted[i].begin();
            phg.doParallelForAllNodes([&](HypernodeID hn){
              if(phg.partID(hn) == from){
                queue.checkAll(hn, i);
              }
              else{
                queue.check(hn, i, from);
              }
            });
          }
        }
        else{
          for(int j = indices[i][from]; j < nodes_sorted[0].size(); j++){
              HypernodeID hn = phg.partWeight(nodes_sorted[i][j].second);
              if(phg.partID(hn) == from){
                queue.checkAll(hn, i);
              }
              else{
                queue.check(hn, i, from);
              }
            }
          while(indices[i][from] < nodes_sorted[0].size() && nodes_sorted[i][indices[i][from]] <= std::abs(_context.partition.max_part_weights[from][i] - phg.partWeight(from).weights[i])){
              indices[i][from]++;
            }
        }
        if(phg.partWeight(to).weights[i] <= _context.partition.max_part_weights[to][i]){
          for(HypernodeID j = nodes_sorted[0].size() - 1; j >= indices[i][from]; j--){
              HypernodeID hn = phg.partID(nodes_sorted[i][j].second);
              if(phg.partID(hn) == to){
                queue.checkAll(hn, i);
              }
            }
            while(indices[i][to] > 0 && nodes_sorted[i][indices[i][to] - 1] > std::abs(_context.partition.max_part_weights[to][i] - phg.partWeight(to).weights[i])){
              indices[i][to]--;
              HypernodeID hn = phg.partID(nodes_sorted[i][indices[i][to]].second);
              if(phg.partID(hn) == from){
                queue.checkAll(hn, i);
              }
              else{}
            }
        }*/
      }

      /*phg.doParallelForAllNodes([&](HypernodeID hn){
        if(phg.partID(hn) == move.from || phg.partID(hn) == move.to){
          for(PartitionID p = 0; p < phg.k(); p++){
            queue.changeBalance({hn, {p, balance_gain(phg, hn, phg.partID(hn), p)}});
          }
        }
        else{
          if(phg.partID(hn) != move.to){
            queue.changeBalance({hn, {move.to, balance_gain(phg, hn, phg.partID(hn), move.to)}});
          }
          if(phg.partID(hn) != move.from){
            queue.changeBalance({hn, {move.from, balance_gain(phg, hn, phg.partID(hn), move.from)}});
          }
        }
      });*/
      /*tbb::parallel_for(UL(0), changed_nodes.size(), [&](HypernodeID hn){
        queue.update(changed_nodes[hn]);
      });*/
      for(size_t i = 0; i < changed_nodes.size(); i++){
        queue.update(changed_nodes[i]);
      }
      queue.checkSizes(phg.k());
    }

      // TODO: rebalancing logic goes here
    
    

    if (moves_by_part != nullptr) {
      moves_by_part->resize(_context.partition.k);
      for (auto& direction : *moves_by_part) direction.clear();
      // TODO: ignore for now, implementation necessary to support unconstrained refinement (compare advanced_rebalancer)
    } else if (moves_linear != nullptr) {
      moves_linear->clear();
      // TODO: ignore for now, implementation necessary to support unconstrained refinement (compare advanced_rebalancer)
    }

    // Update metrics statistics
    Gain delta = _gain.delta();  // note: only correct if all moves were performed with objective_delta as defined above
     HEAVY_REFINEMENT_ASSERT(best_metrics.quality + delta == metrics::quality(phg, _context),
      V(best_metrics.quality) << V(delta) << V(metrics::quality(phg, _context)));
    best_metrics.quality += delta;
    best_metrics.imbalance = metrics::imbalance(phg, _context);

    /*bool improvement = delta < 0;
    return improvement;*/
    return false;
  }
  }


  template <typename GraphAndGainTypes>
  void MDRebalancer<GraphAndGainTypes>::initializeImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph) {
    auto& phg = utils::cast<PartitionedHypergraph>(hypergraph);

    if (!_gain_cache.isInitialized()) {
      _gain_cache.initializeGainCache(phg);
    }
  }

  template <typename GraphAndGainTypes>
  bool MDRebalancer<GraphAndGainTypes>::refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                  const vec<HypernodeID>& , Metrics& best_metrics, double) {
    return refineInternal(hypergraph, nullptr, nullptr, best_metrics);
  }

  template <typename GraphAndGainTypes>
  bool MDRebalancer<GraphAndGainTypes>::refineAndOutputMovesImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                                    const vec<HypernodeID>& ,
                                                                    vec<vec<Move>>& moves_by_part,
                                                                    Metrics& best_metrics,
                                                                    const double) {
    return refineInternal(hypergraph, &moves_by_part, nullptr, best_metrics);
  }

  template <typename GraphAndGainTypes>
  bool MDRebalancer<GraphAndGainTypes>::refineAndOutputMovesLinearImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                                          const vec<HypernodeID>& ,
                                                                          vec<Move>& moves,
                                                                          Metrics& best_metrics,
                                                                          const double) {
    return refineInternal(hypergraph, nullptr, &moves, best_metrics);
  }

  // explicitly instantiate so the compiler can generate them when compiling this cpp file
  namespace {
  #define MD_REBALANCER(X) MDRebalancer<X>
  }

  // explicitly instantiate so the compiler can generate them when compiling this cpp file
  INSTANTIATE_CLASS_WITH_VALID_TRAITS(MD_REBALANCER)
}

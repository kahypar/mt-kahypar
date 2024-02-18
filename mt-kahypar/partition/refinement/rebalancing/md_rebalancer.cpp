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

#include <chrono>

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
                                                      
    auto start = std::chrono::high_resolution_clock::now();                                                  
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
    resizeDataStructuresForCurrentK();
    // This function is passed as lambda to the changeNodePart function and used
    // to calculate the "real" delta of a move (in terms of the used objective function).
    auto objective_delta = [&](const SynchronizedEdgeUpdate& sync_update) {
      _gain.computeDeltaForHyperedge(sync_update);
    };

    // If partition is imbalanced, rebalancer is activated
    if ( !metrics::isBalanced(phg, _context) ) {
      DBG << "Starting multi-dimensional rebalancer";  // only printed if debug=true in header
      _gain.reset();

    auto horizontal_balance_gain = [&](const PartitionedHypergraph& phg, HypernodeID node, PartitionID from, PartitionID to){
      double gain = 0.0;
      if(from == to){
        return 0.0;
      }
      for(int i = 0; i < dimension; i++){
        int32_t to_excess = std::max(0, std::min(phg.nodeWeight(node).weights[i], phg.partWeight(to).weights[i] + 
        phg.nodeWeight(node).weights[i] - _context.partition.max_part_weights[to].weights[i]));
        int32_t from_excess = std::max(0, std::min(phg.nodeWeight(node).weights[i], phg.partWeight(from).weights[i] - _context.partition.max_part_weights[from].weights[i]));
        if(to_excess != 0){
          gain += static_cast<double>(to_excess) * _context.partition.max_part_weights_inv[to][i]; 
        }
        if(from_excess != 0){
          gain -= static_cast<double>(from_excess) * _context.partition.max_part_weights_inv[from][i]; 
        }
      }
      
      return gain;
    };
    auto vertical_imbalance = [&](HypernodeWeight partition_weight, PartitionID p){
      std::vector<double> weight_fractions;
      double sum = 0.0;
      for(int i = 0; i < dimension; i++){
        double fraction = static_cast<double>(partition_weight.weights[i]) * _context.partition.max_part_weights_inv[p][i];
        weight_fractions.push_back(fraction);
        sum += fraction;
      }
      double optimal_fraction = sum / static_cast<double>(dimension);
      double result = 0.0;
      for(int i = 0; i < dimension; i++){
        result += std::abs(weight_fractions[i] - optimal_fraction);
      }
      return result;
    }; 
    auto vertical_balance_gain = [&](const PartitionedHypergraph& phg, HypernodeID node, PartitionID from, PartitionID to){
      return vertical_imbalance(phg.partWeight(from) - phg.nodeWeight(node), from) - vertical_imbalance(phg.partWeight(from), from) 
      + vertical_imbalance(phg.partWeight(to) + phg.nodeWeight(node), to) - vertical_imbalance(phg.partWeight(to), to);
    };

    auto balance_gain = [&](const PartitionedHypergraph& phg, HypernodeID node, PartitionID from, PartitionID to, bool use_horizontal){
      if(use_horizontal){
        return horizontal_balance_gain(phg, node, from, to);
      }
      return vertical_balance_gain(phg, node, from, to);
    };


    auto balance_gain_unweighed = [&](const PartitionedHypergraph& phg, HypernodeID node, PartitionID from, PartitionID to){
      double gain = 0.0;
      if(from == to){
        return 0.0;
      }
      for(int i = 0; i < dimension; i++){
        int32_t to_excess = std::max(0, std::min(phg.nodeWeight(node).weights[i], phg.partWeight(to).weights[i] + 
        phg.nodeWeight(node).weights[i] - _context.partition.max_part_weights[to].weights[i]));
        int32_t from_excess = std::max(0, std::min(phg.nodeWeight(node).weights[i], phg.partWeight(from).weights[i] - _context.partition.max_part_weights[from].weights[i]));
        if(to_excess != 0){
          gain += to_excess; 
        }
        if(from_excess != 0){
          gain -= from_excess; 
        }
      }
      
      return gain;
    }; 


    auto weighed_imbalance = [&](){
      double ib = 0;
      for(PartitionID k = 0; k < phg.k(); k++){
        for(int i = 0; i < dimension; i++){
          ib += std::max(0, phg.partWeight(k).weights[i] - _context.partition.max_part_weights[k].weights[i]) * _context.partition.max_part_weights_inv[k][i];
        }
      }
      return ib;
      
    };

    auto imbalance = [&](){
      int ib = 0;
      for(PartitionID k = 0; k < phg.k(); k++){
        for(int i = 0; i < dimension; i++){
          ib += std::max(0, phg.partWeight(k).weights[i] - _context.partition.max_part_weights[k].weights[i]);
        }
      }
      return ib;
    };
       
    std::vector<Gain> qualities;
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
    std::vector<double> weighed_imbalances;
    std::vector<int> imbalances;
    std::vector<std::vector<HypernodeID>> indices;
    nodes_sorted.resize(mt_kahypar::dimension);
    indices.resize(mt_kahypar::dimension);
    for(int i = 0; i < dimension; i++){
      indices[i].resize(phg.k());
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
      if(nodes_sorted[dimension][middle].first <= std::abs(nw.weights[dimension])){
        left = middle;
      }
      else{
        right = middle;
      }
      }
      return right;
    };
    tbb::parallel_for(UL(0), mt_kahypar::dimension, [&](const size_t i){
      tbb::parallel_for(UL(0), (size_t)phg.k(), [&](const size_t j){
        indices[i][j] = getBorder(i, exceed[j]);
      });     
    });
    
    Gain quality = metrics::quality(phg, _context);
    int num_loops = 0;
    int num_moves = 0;
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
    std::vector<std::vector<std::pair<int32_t, int32_t>>> min_weights;
    std::vector<std::vector<std::pair<int32_t, int32_t>>> max_weights;
    min_weights.resize(phg.k());
    max_weights.resize(phg.k());
    for(PartitionID p = 0; p < phg.k(); p++){
      for(int i = 0; i < dimension; i++){
        min_weights[p].push_back({phg.partWeight(p).weights[i], (int32_t)0});
        max_weights[p].push_back({phg.partWeight(p).weights[i], (int32_t)0});
      }      
    }
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

    int counter = 0;
    //just relevant for fallback case
    double vertical_ib;
    bool horizontal_balance_used = true;
    while(imbalanced != 0 && !queue.isEmpty()){
      num_loops++;
      std::pair<HypernodeID, std::pair<PartitionID, Move_internal>> max_move = queue.getMax();
      ASSERT(max_move.first < phg.initialNumNodes());
      ASSERT(phg.partID(max_move.first) != -1);
      ASSERT(max_move.second.first != -1);
      HypernodeID node = max_move.first;
      PartitionID from = phg.partID(node);
      PartitionID to = max_move.second.first;
      double balance = balance_gain(phg, node, from, to, horizontal_balance_used);
      if(max_move.second.second.balance > 0.00000001){
        std::cout << max_move.second.second.balance << "\n";
      }
      if(balance > max_move.second.second.balance + 0.00000000001){
        queue.changeBalance({node, {to, balance}});
        queue.update(node);
      }
      else{
        num_moves++;
      queue.deleteMax();
      counter++;
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
      if(horizontal_balance_used){
        queue.lock(node);
      }
      tbb::concurrent_vector<HypernodeID> changed_nodes = _gain.getChangedMoves(phg, {move.from, move.to, move.node, 0}, &queue);
      for(PartitionID p = 0; p < phg.k(); p++){
        queue.addToGain({node, {p, -gains.gain}});
        queue.changeBalance({node, {p, balance_gain(phg, node, to, p, horizontal_balance_used)}});
      }
      changed_nodes.push_back(node);

      for(int i = 0; i < dimension; i++){
        if(exceed[from].weights[i] < 0 && (exceed_from_old.weights[i] > 0 || -exceed_from_old.weights[i] < nodes_sorted[i][nodes_sorted[i].size() - 1].first)){
          min_weights[from][i] = {phg.partWeight(from).weights[i], counter};
          HypernodeID border = exceed_from_old.weights[i] < 0 ? getBorder(i, exceed_from_old) : 0;
          for(HypernodeID hn = border; hn < nodes_sorted[i].size(); hn++){
            if(phg.partID(nodes_sorted[i][hn].second) != from){
              queue.changeBalance({nodes_sorted[i][hn].second, {from, balance_gain(phg, nodes_sorted[i][hn].second, 
                phg.partID(nodes_sorted[i][hn].second), from, horizontal_balance_used)}});
              changed_nodes.push_back(nodes_sorted[i][hn].second);
            }            
          }         
        }
        if(exceed[to].weights[i] > 0 && (exceed_to_old.weights[i] < 0 ||exceed_to_old.weights[i] < nodes_sorted[i][nodes_sorted[i].size() - 1].first)){
          max_weights[to][i] = {phg.partWeight(to).weights[i], counter};
          HypernodeID border = exceed_to_old.weights[i] > 0 ? getBorder(i, exceed_to_old) : 0;
          for(HypernodeID hn = border; hn < nodes_sorted[i].size(); hn++){
            if(phg.partID(nodes_sorted[i][hn].second) == to){
              for(PartitionID p = 0; p < phg.k(); p++){
                if(p != to){
                  queue.changeBalance({nodes_sorted[i][hn].second, {p, balance_gain(phg, nodes_sorted[i][hn].second, 
                  phg.partID(nodes_sorted[i][hn].second), p, horizontal_balance_used)}});
                  changed_nodes.push_back(nodes_sorted[i][hn].second);
                }                
              }              
            }            
          }  
        }
      }
      for(size_t i = 0; i < changed_nodes.size(); i++){
        queue.update(changed_nodes[i]);
      }
        weighed_imbalances.push_back(weighed_imbalance());
        imbalances.push_back(imbalance());
        qualities.push_back(best_metrics.quality + _gain.delta());
      queue.checkSizes(phg.k());
      if(!horizontal_balance_used){
        double current_vertical_ib = 0.0;
        for(PartitionID p = 0; p < phg.k(); p++){
          current_vertical_ib += vertical_imbalance(phg.partWeight(p), p); 
        }
        if(current_vertical_ib <= 0.7 * vertical_ib){
          horizontal_balance_used = true;
          for(HypernodeID hn = 0; hn < phg.initialNumNodes(); hn++){
            for(PartitionID p = 0; p < phg.k(); p++){
              queue.changeBalance({hn, {p, balance_gain(phg, hn, phg.partID(hn), p, horizontal_balance_used)}});
            }
          }
          for(HypernodeID hn = 0; hn < phg.initialNumNodes(); hn++){
            queue.update(hn);
          }
        }
      }
      }
      
      if(imbalanced != 0 && queue.isEmpty()){
        std::cout << "\n\n\n\n\n\n\nfallback activated";
        horizontal_balance_used = false;
        for(HypernodeID hn = 0; hn < phg.initialNumNodes(); hn++){
          for(PartitionID p = 0; p < phg.k(); p++){
            queue.changeBalance({hn, {p, balance_gain(phg, hn, phg.partID(hn), p, horizontal_balance_used)}});
          }
        }
        for(HypernodeID hn = 0; hn < phg.initialNumNodes(); hn++){
          queue.update(hn);
        }
        vertical_ib = 0.0;
        for(PartitionID p = 0; p < phg.k(); p++){
          vertical_ib += vertical_imbalance(phg.partWeight(p), p); 
        }

        /*PartitionID max_partition = 0;
        int max_dimension = 0;
        double max_exceed = 0.0;
        for(int dim = 0; dim < dimension; dim++){
          for(PartitionID p = 0; p < phg.k(); p++){
            double exceed = static_cast<double>(phg.partWeight(p).weights[dim]) / static_cast<double>(_context.partition.max_part_weights[p][dim]);
            if(exceed - 1.0 > max_exceed){
              max_exceed = exceed - 1.0;
              max_dimension = dim;
              max_partition = p;
            }
          }
        }
        std::vector<std::vector<HypernodeID>> exchange_nodes;
        exchange_nodes*/
        
      }
    }
    std::cout << queue.isEmpty() << " " << imbalanced << "\n\n\n" << (imbalanced != 0 && queue.isEmpty()) << "\n";
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
    ASSERT(best_metrics.quality + delta == metrics::quality(phg, _context));
     HEAVY_REFINEMENT_ASSERT(best_metrics.quality + delta == metrics::quality(phg, _context),
      V(best_metrics.quality) << V(delta) << V(metrics::quality(phg, _context)));
    best_metrics.quality += delta;
    best_metrics.imbalance = metrics::imbalance(phg, _context);
    auto end = std::chrono::high_resolution_clock::now();  
    std::chrono::duration<double> elapsed_seconds = end - start;
    auto x = std::chrono::duration_cast<std::chrono::microseconds>(elapsed_seconds);

    /*std::cout << "\nRebalancing Round: " << phg.initialNumNodes() << " " << x.count() << " " << num_loops << " " << num_moves << " " << queue.top_moves.swap_down_counter << " " <<
      queue.top_moves.swap_up_counter << " " << quality << " " << metrics::quality(phg, _context) << "\n";
    std::cout << "Weighed Imbalances:\n";
    for(int i = 0; i < weighed_imbalances.size(); i++){
      std::cout << weighed_imbalances[i] << " " ;
    }
    std::cout << "\n";
    std::cout << "Imbalances:\n";
    for(int i = 0; i < imbalances.size(); i++){
      std::cout << imbalances[i] << " " ;
    }
    std::cout << "\n";
    std::cout << "Gains:\n";
    for(int i = 0; i < qualities.size(); i++){
      std::cout << qualities[i] << " " ;
    }
    std::cout << "\n";
    std::cout << imbalanced << "\n";
    for(HypernodeID hn : phg.nodes()){
      for(PartitionID p = 0; p < phg.k(); p++){
        if(balance_gain(phg, hn, phg.partID(hn), p, true) < 0.0){
          std::cout << "\n\n\n\nerror" << balance_gain(phg, hn, phg.partID(hn), p, true);
        }
        if(phg.partWeight(phg.partID(hn)) > _context.partition.max_part_weights[phg.partID(hn)]){
          std::cout << phg.nodeWeight(hn).weights[0] << " " << phg.nodeWeight(hn).weights[1] << " " <<balance_gain(phg, hn, phg.partID(hn), p, true) << "\n";
        }
      }
    }*/
    
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

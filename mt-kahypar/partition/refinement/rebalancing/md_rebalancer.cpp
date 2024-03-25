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
#include <list>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/cast.h"

#include "mt-kahypar/parallel/parallel_prefix_sum.h"

#include <queue>

namespace mt_kahypar{


  struct interval{
    HypernodeID lower_index;
    HypernodeID upper_index;
    double lower_search_value;
    double upper_search_value;
  };

  struct pointer_list{
    std::vector<std::pair<int,HypernodeID>> next_per_level;
  };


  struct hypernodes_ordered{
    std::vector<std::vector<double>> nodes_normalized;
    std::vector<std::pair<double,HypernodeID>> nodes_normalized_total_weight;
    std::vector<int> imbalanced_dimensions;
    std::vector<std::vector<std::pair<HypernodeID, double>>> nodes_prio;
    std::vector<std::vector<pointer_list>> pointerList;
    std::vector<HypernodeID> current_indices;


    //dimensions has to be sorted
    void initialize(std::vector<std::vector<double>> nn, std::vector<std::pair<double,HypernodeID>> nntw, std::vector<bool> could_be_imbalanced){
      nodes_normalized = nn;
      nodes_normalized_total_weight = nntw;
      for(int i = 0; i < could_be_imbalanced.size(); i++){
        if(could_be_imbalanced[i]){
          imbalanced_dimensions.push_back(i);
        }
      }
      nodes_prio.resize(std::pow(2, imbalanced_dimensions.size()));
      current_indices.resize(std::pow(2, imbalanced_dimensions.size()), 0);
      pointerList.resize(std::pow(2, imbalanced_dimensions.size()));
    }

    void initialize(std::vector<std::vector<double>> nn, std::vector<std::pair<double,HypernodeID>> nntw){
      nodes_normalized = nn;
      nodes_normalized_total_weight = nntw;
      for(int i = 0; i < dimension; i++){
        imbalanced_dimensions.push_back(i);
      }
      nodes_prio.resize(std::pow(2, imbalanced_dimensions.size()));
      pointerList.resize(std::pow(2, imbalanced_dimensions.size()));
      current_indices.resize(std::pow(2, imbalanced_dimensions.size()), 0);
    }

    size_t get_index(std::vector<bool> imbalanced){
      size_t index = 0;
      for(int i = 0; i < imbalanced_dimensions.size(); i++){
        if(imbalanced[imbalanced_dimensions[i]]){
          index += std::pow(2,i);
        }
      }
      return index;
    }

    //imbalanced has to be sorted
    std::vector<std::pair<HypernodeID, double>> nodes_priorised(std::vector<bool> imbalanced){
      size_t index = get_index(imbalanced);
      pointerList[index].resize(nodes_normalized.size() + 1);
      std::vector<size_t> imbalanced_list;
      for(int i = 0; i < imbalanced_dimensions.size(); i++){
        if(imbalanced[imbalanced_dimensions[i]]){
          imbalanced_list.push_back(imbalanced_dimensions[i]);
        }
      }
      if(nodes_prio[index].size() != 0){
        return nodes_prio[index];
      }
      nodes_prio[index].resize(nodes_normalized.size());
      for(HypernodeID hn = 0; hn < nodes_normalized.size(); hn++){
        nodes_prio[index][hn] = {hn, 0.0};
        for(size_t i = 0; i < imbalanced_list.size(); i++){
          nodes_prio[index][hn].second += nodes_normalized[hn][imbalanced_list[i]];
        }
        nodes_prio[index][hn].second /= nodes_normalized_total_weight[hn].first;
      }
      std::sort(nodes_prio[index].begin(), nodes_prio[index].end(), [&](std::pair<HypernodeID, double> a, std::pair<HypernodeID, double> b){
        return a.second > b.second;
      });
      std::list<std::pair<int,HypernodeID>> current_last;
      for(HypernodeID hn = nodes_prio[index].size() - 1; hn-- > 0;){
        HypernodeID node = nodes_prio[index][hn].first;
        int level_number = 0;
        double ignored_percentage = 0.1;
        while(std::ceil(ignored_percentage * nodes_normalized.size()) < nodes_normalized_total_weight[node].second + 1){
          level_number++;
          ignored_percentage = (ignored_percentage + 1.0) / 2.0;
        }
        int idx_in_current = 0;

        for(auto x : current_last){
          pointerList[index][hn + 1].next_per_level.push_back(x);
        }       
        
        auto it = current_last.begin();
        while(it != current_last.end() && (*it).first < level_number){
          std::next(it, 1);
        }
        current_last.erase(it, current_last.end());
        
        current_last.push_back({level_number, hn + 1});
              
      }
      
      
      for(auto x : current_last){
        pointerList[index][0].next_per_level.push_back(x);
      }






      /*std::vector<pointer_list> &p_l = pointerList[index];
      double frac = 0.5;
      HypernodeID level_index;
      int level_number = 0;
      std::vector<bool> nodeInserted(nodes_prio[index].size(), false);
      do{
        level_index = std::ceil(frac * nodes_normalized.size());
        HypernodeID lastIndex = 0;
        for(HypernodeID hn = 0; hn < nodes_prio[index].size(); hn++){
          if(nodes_normalized_total_weight[nodes_prio[index][hn].first].second <= level_index){
            if(!nodeInserted[hn]){
              p_l[lastIndex].next_per_level.push_back({level_number, hn});
            }
            lastIndex = hn + 1;
            nodeInserted[hn] = true;
          }
        }
        frac = (1.0 + frac) / 2.0;
        level_number++;
      }while(level_index != nodes_normalized.size());
      p_l[p_l.size() - 1].next_per_level.push_back({0, -1});*/
      
      return nodes_prio[index];
    }

    std::pair<HypernodeID,int> getNext(std::vector<bool> imbalanced, int level){
      nodes_priorised(imbalanced);
      size_t index = get_index(imbalanced);
      pointer_list current = pointerList[index][current_indices[index]];
      if(current.next_per_level.size() == 0){
        return {0, -1};
      }
      int start = 0;
      int end = current.next_per_level.size();
      while(end != start + 1 && current.next_per_level[start + 1].first < level){
        int middle = (start + end) / 2;
        if(current.next_per_level[middle].first <= level){
          start = middle;
        }
        else{
          end = middle;
        }
      }
      current_indices[index] = current.next_per_level[start].second;
      int lowest_above = start < current.next_per_level.size() - 1 ? current.next_per_level[start + 1].first : current.next_per_level[start].first;
      return {nodes_prio[index][current.next_per_level[start].second - 1].first, lowest_above};
    }


    void reset(){
      for(HypernodeID hn = 0; hn < current_indices.size(); hn++){
        current_indices[hn] = 0;
      }
    }


  };

  /*std::vector<double> relative_weights(HypernodeWeight nw, HypernodeWeight max){
    std::vector<double> res;
    double sum = 0.0;
    for(int i = 0; i < dimension; i++){
      double fraction = static_cast<double>(nw.weights[i]) / static_cast<double>(max.weights[i]);
      res.push_back(fraction);
      sum += fraction;
    }
    for(int i = 0; i < dimension; i++){
      res[i] /= sum;
    }
    return res;
  }*/


  double vertical_imbalance(std::vector<double> weights){
    double total = 0.0;
    for(int i = 0; i < weights.size(); i++){
      total += weights[i];
    }
    total /= weights.size();
    double loss = 0.0;
    for(int i = 0; i < weights.size(); i++){
      loss += (weights[i] - total) * (weights[i] - total);
    }
    return loss;
  }

  std::vector<double> minus(std::vector<double> a, std::vector<double> b){
    std::vector<double> res;
    for(int i = 0; i < a.size(); i++){
      res.push_back(a[i] - b[i]);
    }
    return res;
  }

  std::pair<bool, std::vector<std::vector<std::pair<PartitionID, HypernodeID>>>> bin_packing(std::vector<std::vector<double>> limits, std::vector<std::pair<double, std::pair<PartitionID, HypernodeID>>> nodes, std::vector<std::vector<std::vector<double>>> weights){
    std::vector<std::vector<std::pair<PartitionID, HypernodeID>>> packing(limits.size());
    std::cout << "binpack\n";
    for(HypernodeID hn = 0; hn < nodes.size(); hn++){
      double max_bg = std::numeric_limits<double>::max();
      PartitionID max_p = -1;
      for(PartitionID p = 0; p < limits.size(); p++){
        std::vector<double> new_limit = minus(limits[p], weights[nodes[hn].second.first][nodes[hn].second.second]);
        bool too_big = false;
        for(int d = 0; d < dimension; d++){
          if(new_limit[d] <= 0.0){
            too_big = true;
          }
        }
        if(too_big){
          continue;
        }
        double bg = vertical_imbalance(new_limit) - vertical_imbalance(limits[p]);
        if(bg < max_bg){
          max_p = p;
          max_bg = bg;
        }
      }
      if(max_p == -1){
        return {false, packing};
      }
      limits[max_p] = minus(limits[max_p], weights[nodes[hn].second.first][nodes[hn].second.second]);
      packing[max_p].push_back(nodes[hn].second);
    }
    return {true, packing};
  }




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
      Gain isolated_block_gain = 0;
      std::vector<Move_md> moves = _gain.allGains(phg, hn, isolated_block_gain);
      for(int i = 0; i < moves.size(); i++){
        ASSERT(phg.partID(hn) != moves[i].to);
        ASSERT(phg.partID(hn) != -1);
        queue.insert_without_updating({hn, {moves[i].to, {moves[i].gain_and_balance, moves[i].gain, moves[i].balance}}});
      }
    });
    queue.check();
    PartitionID imbalanced = 0;
    for(PartitionID p = 0; p < phg.k(); p++){
      if(phg.partWeight(p) > _context.partition.max_part_weights[p]){
        imbalanced++;
      }
    }

    int counter = 0;
    //just relevant for fallback case
    double vertical_ib;
    std::vector<Fallback_MoveQueue> moves_per_partition;
    moves_per_partition.resize(phg.k());
    bool horizontal_balance_used = true;
    std::vector<HypernodeID> min_affected_node;
    for(int d = 0; d < dimension; d++){
      min_affected_node.push_back(phg.initialNumNodes());
    } 
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
      if(!horizontal_balance_used){
        queue.lock(node);
        queue.update(node);
      }
      tbb::concurrent_vector<HypernodeID> changed_nodes = _gain.getChangedMoves(phg, {move.from, move.to, move.node, 0}, &queue);
      for(PartitionID p = 0; p < phg.k(); p++){
        queue.addToGain({node, {p, -gains.gain}});
        queue.changeBalance({node, {p, balance_gain(phg, node, to, p, horizontal_balance_used)}});
      }
      changed_nodes.push_back(node);

      for(int i = 0; i < dimension; i++){
        if(exceed[from].weights[i] < 0 && (exceed_from_old.weights[i] > 0 || -exceed_from_old.weights[i] < nodes_sorted[i][nodes_sorted[i].size() - 1].first)){
          HypernodeID border = exceed_from_old.weights[i] < 0 ? getBorder(i, exceed_from_old) : 0;
          min_affected_node[i] = std::min(min_affected_node[i], border);
        }
        if(exceed[to].weights[i] > 0 && (exceed_to_old.weights[i] < 0 || exceed_to_old.weights[i] < nodes_sorted[i][nodes_sorted[i].size() - 1].first)){
          HypernodeID border = exceed_to_old.weights[i] > 0 ? getBorder(i, exceed_to_old) : 0;
          min_affected_node[i] = std::min(min_affected_node[i], border);
        }
        if(counter % 50 == 0){
          for(HypernodeID hn = min_affected_node[i]; hn < nodes_sorted[i].size(); hn++){
            for(PartitionID p = 0; p < phg.k(); p++){
              if(p != phg.partID(hn)){
                queue.changeBalance({nodes_sorted[i][hn].second, {p, balance_gain(phg, nodes_sorted[i][hn].second, 
                  phg.partID(nodes_sorted[i][hn].second), p, horizontal_balance_used)}});
              }
            }
            changed_nodes.push_back(hn);
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
        if(current_vertical_ib <= 0.9 * vertical_ib){
          horizontal_balance_used = true;
          for(PartitionID p = 0; p < phg.k(); p++){
            if(phg.partWeight(p) > _context.partition.max_part_weights[p]){
              std::cout << p << " ib\n";
            }
            for(int j = 0; j < dimension; j++){
            std::cout << phg.partWeight(p).weights[j] << " ";
          }
          std::cout << "\n";
          std::cout << vertical_imbalance(phg.partWeight(p), p) << "\n"; 
          }
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
        std::cout << "fb activated" << "\n";
        /*std::vector<HypernodeID> L;
        std::vector<HypernodeID> V_N;
        std::vector<std::pair<HypernodeID, Gain>> prio_for_S;
        std::vector<std::pair<HypernodeID, Gain>> node_imbalance;
        for(HypernodeID hn : phg.nodes()){
          if(phg.nodeWeight(hn).max() <= _context.partition.epsilon[0]){
            L.push_back(hn);
          }
          else{
            V_N.push_back(hn);
          }
        }
        for(HypernodeID hn : V_N){
          Gain isolated_block_gain = 0;
          _gain.allGains(phg, hn, isolated_block_gain);
          prio_for_S.push_back({hn, isolated_block_gain / (phg.nodeWeight(hn).max() - phg.nodeWeight().min() + 0.1)});
          node_imbalance.push_back({hn, phg.nodeWeight(hn).max() - phg.nodeWeight().min()});
        }
        std::sort(prio_for_S.begin(), prio_for_S.end(), [&](std::pair<Hypernode,Gain> a, std::pair<Hypernode,Gain> b){
          return a.second <= b.second;
        });
        std::sort(node_imbalance.begin(), node_imbalance.end(), [&](std::pair<Hypernode,Gain> a, std::pair<Hypernode,Gain> b){
          return a.second <= b.second;
        });
        std::vector<HypernodeWeight> current_partition_weights;
        for(PartitionID p = 0; p < phg.k(); p++){
          current_partition_weights.push_back(phg.partWeight(p));
        }
        auto remove_nodes 
        HypernodeID start = 0;
        HypernodeID end = 1;
        bool balanced = false;*/

      auto smaller_equal = [](std::vector<double> a, std::vector<double> b){
        for(size_t i = 0; i < a.size(); i++){
          if(a[i] > b[i]){
            return false;
          }
        }
        return true;
      };





      std::vector<HypernodeID> L;
      const double L_threshold = 0.0;

      std::vector<std::vector<std::vector<double>>> normalized_weights(phg.k());
      //total Weight of each node (organized per partition)
      std::vector<std::vector<std::pair<double,HypernodeID>>> totalWeight(phg.k());
      std::vector<std::vector<std::pair<double,HypernodeID>>> totalWeight_ordered(phg.k());
      std::vector<std::vector<double>> normalized_partition_weights(phg.k());
      std::vector<std::vector<double>> normalized_perfect_weights(phg.k());
      std::vector<bool> is_imbalanced(phg.k(), false);
      std::vector<std::pair<PartitionID, bool>> imbalances_per_partition(phg.k());
      std::vector<struct hypernodes_ordered> nodes_ordered(phg.k());
      std::vector<HypernodeID> map_id_to_index(phg.initialNumNodes());
      std::vector<std::vector<HypernodeID>> map_index_to_id(phg.k());
      

      
      for(int p = 0; p < phg.k(); p++){
        normalized_partition_weights[p].resize(dimension);
        normalized_perfect_weights[p].resize(dimension);
        for(int d = 0; d < dimension; d++){
          normalized_partition_weights[p][d] = phg.partWeight(p).weights[d] * _context.partition.max_part_weights_inv[p][d];
          normalized_perfect_weights[p][d] = _context.partition.perfect_balance_part_weights[p].weights[d] * _context.partition.max_part_weights_inv[p][d];
        }
      }

      //initialize normalized_weights, totalWeight and normalized_partition_weights
      for(HypernodeID hn : phg.nodes()){
        PartitionID part = phg.partID(hn);
        double total_weight = 0.0;
        std::vector<double> weight;
        double max_weight = 0.0;
        for(int d = 0; d < dimension; d++){
          double w = phg.nodeWeight(hn).weights[d] * _context.partition.max_part_weights_inv[phg.partID(hn)][d];
          total_weight += w;
          weight.push_back(w);
          max_weight = std::max(max_weight, w);         
        }
        if(max_weight > L_threshold){
          map_id_to_index[hn] = normalized_weights[part].size();
          map_index_to_id[part].push_back(hn);
          normalized_weights[part].push_back(weight);
          totalWeight_ordered[part].push_back({total_weight, map_id_to_index[hn]});
        }
        else{
          L.push_back(hn);
        }
        
      }


      //initialize imbalances_per_partition, is_imbalanced
      /*for(PartitionID p = 0; p < phg.k(); p++){
        imbalances_per_partition[p].resize(dimension, false);
        for(int j = 0; j < dimension; j++){
          if(phg.partWeight(p).weights[j] > _context.partition.max_part_weights[p].weights[j]){
            is_imbalanced[p] = true;
            imbalances_per_partition[p][j] = true;
          }
        }
      }*/


      for(PartitionID p = 0; p < phg.k(); p++){
        std::sort(totalWeight_ordered[p].begin(), totalWeight_ordered[p].end());
        totalWeight[p].resize(totalWeight_ordered[p].size());
        for(HypernodeID hn = 0; hn < totalWeight_ordered[p].size(); hn++){
          totalWeight[p][(totalWeight_ordered[p][hn]).second] = {(totalWeight_ordered[p][hn]).first, hn};
        }
      }


      //initialize nodes_ordered
      for(PartitionID p = 0; p < phg.k(); p++){
        nodes_ordered[p].initialize(normalized_weights[p], totalWeight[p]);
      }
      



      std::vector<std::pair<PartitionID, double>> extraction_order;
      std::vector<std::vector<double>> max_remaining_weight(phg.k());
      double factor_lowerBound = 0.0;
      /*double factor_upperBound = 1.0 + *std::min_element(_context.partition.epsilon.begin(), _context.partition.epsilon.end());*/
      double factor_maxbound = 1.0 + _context.partition.epsilon[0];
      double factor_upperBound = factor_maxbound;
      double factor = 1.0;
      double search_threshold = 0.1;
      double fraction_nodes = 0.5;
      std::vector<int> current_level_per_partition(phg.k(), 0);
      std::vector<std::vector<std::pair<PartitionID, HypernodeID>>> current_result;
      for(PartitionID p = 0; p < phg.k(); p++){
        double excess = 0.0;
        for(int d = 0; d < dimension; d++){
          excess += std::max(0.0, normalized_partition_weights[p][d] - normalized_perfect_weights[p][d]);
        }
        extraction_order.push_back({p, excess});
      }
      std::sort(extraction_order.begin(), extraction_order.end(), [&](std::pair<PartitionID, double> a, std::pair<PartitionID, double> b){
        return a.second > b.second;
      });
      std::cout << "marker3\n";
      while(true){
        for(PartitionID p = 0; p < phg.k(); p++){
          max_remaining_weight[p].resize(dimension);
          for(int d = 0; d < dimension; d++){
            max_remaining_weight[p][d] = normalized_perfect_weights[p][d] * factor;
          }
        }
        bool success = false;
        std::vector<std::vector<double>> virtual_weight(phg.k());
        std::vector<std::vector<bool>> extracted(phg.k());
        for(int i = 0; i < extraction_order.size(); i++){
          PartitionID p = extraction_order[i].first;
          for(int d = 0; d < dimension; d++){
            virtual_weight[p].push_back(normalized_partition_weights[p][d]);
          }
          if(smaller_equal(normalized_partition_weights[p], max_remaining_weight[p])){
            continue;
          }
          bool improvement = true;
          int level = current_level_per_partition[p] - 1;
          double level_quality = std::numeric_limits<double>::max();
          while(improvement){
            std::vector<bool> imbalanced(dimension);
            std::vector<bool> used(normalized_weights[p].size(), false);
            int highest_would_use = 0;
            std::vector<double> this_run_weight;
            for(int d = 0; d < dimension; d++){
              this_run_weight.push_back(normalized_partition_weights[p][d]);
            }
            bool possible_to_balance = false;
            while(!smaller_equal(this_run_weight, max_remaining_weight[p])){
              for(int d = 0; d < dimension; d++){
                imbalanced[d] = this_run_weight[d] > max_remaining_weight[p][d];
              }
              std::pair<HypernodeID,int> next;
              do{
                next = nodes_ordered[p].getNext(imbalanced, level + 1);
              }while(next.second != -1 && used[next.first]);
              if(next.second == -1){
                std::cout <<"minusone\n";
                break;
              }
              used[next.first] = true;
              highest_would_use = std::max(highest_would_use, next.second);
              for(int d = 0; d < dimension; d++){
                this_run_weight[d] -= normalized_weights[p][next.first][d];
              }
              possible_to_balance = true;
            }
            
            nodes_ordered[p].reset();
            if(!possible_to_balance){
              std::cout <<"np\n";
              level++;
              continue;
            }
            if(smaller_equal(this_run_weight, max_remaining_weight[p])){
              double quality = 0.0;
              for(int d = 0; d < dimension; d++){
                quality += max_remaining_weight[p][d] - this_run_weight[d];
              }
              if(quality >= level_quality){
                improvement = false;
                level--;
              }
              else{
                level_quality = quality;
                extracted[p] = used;
                virtual_weight[p] = this_run_weight;
                if(highest_would_use <= level + 1){
                  improvement = false;
                }
              }
            }
            
            
            level++;
          }
          current_level_per_partition[p] = level;
        }
        std::cout << "marker4\n";
        std::vector<std::pair<double, std::pair<PartitionID,HypernodeID>>> bp_nodes;
        for(PartitionID p = 0; p < phg.k(); p++){
          for(HypernodeID hn = 0; hn < extracted[p].size(); hn++){
            if(extracted[p][hn]){
              HypernodeID node = map_index_to_id[p][hn];
              double min = std::numeric_limits<double>::max();
              double max = std::numeric_limits<double>::min();
              for(int d = 0; d < dimension; d++){
                min = std::min(min, normalized_weights[p][hn][d]);
                max = std::max(max, normalized_weights[p][hn][d]);
              }

              bp_nodes.push_back({(max - min) / max, {p, hn}});
            }
          }
        }
        
        std::sort(bp_nodes.begin(), bp_nodes.end(), [&](std::pair<double, std::pair<PartitionID,HypernodeID>> a, 
          std::pair<double, std::pair<PartitionID,HypernodeID>> b){
            return a.first > b.first;
          });
          
        std::vector<std::vector<double>> limits(phg.k());
        for(PartitionID p = 0; p < phg.k(); p++){
          std::vector<double> entry(dimension);
          for(int d = 0; d < dimension; d++){
            entry[d] = 1 - virtual_weight[p][d];
          }
          limits[p] = entry;
        }
        
        std::pair<bool, std::vector<std::vector<std::pair<PartitionID, HypernodeID>>>> result = bin_packing(limits, bp_nodes, normalized_weights);


        if(result.first){
          factor_lowerBound = factor;
          current_result = result.second;
        }
        else{
          std::cout << "false\n";
          factor_upperBound = factor;
        }
        if(factor_upperBound - factor_lowerBound < search_threshold){
          break;
        }
        factor = std::max((factor_upperBound + factor_lowerBound) / 2.0, 2.0 * factor_upperBound - factor_maxbound);
      }

      for(PartitionID p = 0; p < phg.k(); p++){
        for(HypernodeID hn = 0; hn < current_result[p].size(); hn++){
          HypernodeID node = map_index_to_id[current_result[p][hn].first][current_result[p][hn].second];
          phg.changeNodePart(node, phg.partID(node), p, objective_delta);
        }
      }








        

        /*double max = 0.0;
        PartitionID max_part = 0;
        int max_dimension = 0;
        for(PartitionID p = 0; p < phg.k(); p++){
          for(int d = 0; d < dimension; d++){
            if((phg.partWeight(p).weights[d] - _context.partition.max_part_weights[p].weights[d]) * _context.partition.max_part_weights_inv[p][d] > max){
              max = (phg.partWeight(p).weights[d] - _context.partition.max_part_weights[p].weights[d]) * _context.partition.max_part_weights_inv[p][d];
              max_dimension = d;
              max_part = p;
            }
          }
        }
        std::vector<std::vector<std::pair<double, HypernodeID>>> nodes_fraction_list;
        nodes_fraction_list.resize(phg.k());
        for(HypernodeID hn = 0; hn < phg.initialNumNodes(); hn++){
          nodes_fraction_list[phg.partID(hn)].push_back({
            static_cast<double>(phg.nodeWeight(hn).weights[max_dimension]) / static_cast<double>(phg.nodeWeight(hn).weights()), hn});
          
        }
        tbb::parallel_for(UL(0), phg.k(), [&](const size_t i){tbb::parallel_sort(nodes_fraction_list[i].begin(), nodes_fraction_list[i].end());});*/





        /*std::cout << "\nfallback activated\n";
        std::vector<int> numbers;
        numbers.resize(phg.k());
        horizontal_balance_used = false;
          for(PartitionID p = 0; p < phg.k(); p++){
            numbers[p] = 0;
            moves_per_partition[p].initialize(phg.initialNumNodes(), phg.k());
          }
        for(HypernodeID hn = 0; hn < phg.initialNumNodes(); hn++){
          for(PartitionID p = 0; p < phg.k(); p++){
            FallBackMove move = {0.0, queue.get({hn, p}).gain, balance_gain(phg, hn, phg.partID(hn), p, horizontal_balance_used)};
            move.recomputeBalance();
            moves_per_partition[phg.partID(hn)].insert({hn, {p, move}});
            numbers[phg.partID(hn)]++;
          }
        }
        while(imbalanced != 0){
          PartitionID max_partition = 0;
          int max_dimension = 0;
          double max_exceed = 0.0;
          for(int dim = 0; dim < dimension; dim++){
            for(PartitionID p = 0; p < phg.k(); p++){
              double exceed = static_cast<double>(phg.partWeight(p).weights[dim]) / static_cast<double>(_context.partition.max_part_weights[p].weights[dim]);
              if(exceed - 1.0 > max_exceed){
                max_exceed = exceed - 1.0;
                max_dimension = dim;
                max_partition = p;
              }
            }
          }
          if(moves_per_partition[max_partition].isEmpty()){
          }
          else{
            std::cout << max_partition << " " << numbers[max_partition] << "\n";
            numbers[max_partition]--;
            std::pair<HypernodeID, std::pair<PartitionID, FallBackMove>> move = moves_per_partition[max_partition].getMax();
            if(balance_gain(phg, move.first, phg.partID(move.first), move.second.first, horizontal_balance_used) > move.second.second.balance + 0.00001){
              moves_per_partition[max_partition].changeBalance({move.first, {move.second.first, 
                balance_gain(phg, move.first, phg.partID(move.first), move.second.first, horizontal_balance_used)}});
              moves_per_partition[max_partition].update(move.first);
              continue;
            }
          HypernodeID node = move.first;
          PartitionID to = move.second.first;
          imbalanced += (phg.partWeight(phg.partID(node)) - phg.nodeWeight(node) > _context.partition.max_part_weights[phg.partID(node)])
        - (phg.partWeight(phg.partID(node)) > _context.partition.max_part_weights[phg.partID(node)])
        - (phg.partWeight(to) > _context.partition.max_part_weights[to])
        + (phg.partWeight(to) + phg.nodeWeight(node) > _context.partition.max_part_weights[to]);
          phg.changeNodePart(node, phg.partID(node), to, objective_delta);
          tbb::concurrent_vector<HypernodeID> changed_nodes = _gain.getChangedMoves(phg, {max_partition, to, node, 0}, &queue, &moves_per_partition);
          for(HypernodeID hn : changed_nodes){
            moves_per_partition[phg.partID(hn)].update(hn);
          }
          }
          std::vector<std::vector<HypernodeID>> nodes_per_permutation;
          nodes_per_permutation.resize(faculty(dimension));
          for(HypernodeID hn : phg.nodes()){
            nodes_per_permutation[getIndex(getPermutation(phg.nodeWeight(hn)))].push_back(hn);
          }
          
        }*/
        
      }
    }
    /*std::cout << queue.isEmpty() << " " << imbalanced << "\n\n\n" << (imbalanced != 0 && queue.isEmpty()) << "\n";*/
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
    /*ASSERT(best_metrics.quality + delta == metrics::quality(phg, _context));*/
     HEAVY_REFINEMENT_ASSERT(best_metrics.quality + delta == metrics::quality(phg, _context),
      V(best_metrics.quality) << V(delta) << V(metrics::quality(phg, _context)));
    best_metrics.quality += delta;
    best_metrics.imbalance = metrics::imbalance(phg, _context);
    auto end = std::chrono::high_resolution_clock::now();  
    std::chrono::duration<double> elapsed_seconds = end - start;
    auto x = std::chrono::duration_cast<std::chrono::microseconds>(elapsed_seconds);
    
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


  /*std::vector<int> permute(std::vector<int> vec, int x, int y){
    std::vector<int> res;
    for(int i = 0; i < vec.size() - 1; i++){
      res.push_back(vec[i]);
    }
    int tmp = res[x];
    res[x] = res[y];
    res[y] = tmp;
    return res;
  }

  std::vector<int> getPermutation(std::vector<double> weight){
    std::vector<int> permutation;
    for(int i = 0; i < weight.size(); i++){
      permutation.push_back(-1);
    }
    for(int i = 0; i < weight.size(); i++){
      int rank = 0;
      for(int j = 0; j < weight.size(); j++){
        if(weight[i] > weight[j]){
          rank++;
        }
      }
      while(permutation[rank] != -1){
        rank++;
      }
      permutation[rank] = i;
    }
    return permutation;
  }

  std::vector<std::pair<double, std::vector<int>>> weightPermutations(std::vector<double> weight){
    struct element{
      double penalty;
      int first_idx;
      int current_idx;
      std::vector<int> permutation;
    };
    std::vector<std::pair<double, std::vector<int>>> result;
    std::priority_queue<element> queue;
    queue.push({0.0, 0, 1, getPermutation(weight)});
    while(!queue.empty()){
      element next = queue.top();
      result.push_back({next.penalty, next.permutation});
      queue.pop();
      if(next.first_idx != 0){
        queue.push({weight[next.permutation[next.first_idx]] - weight[next.permutation[next.first_idx - 1]], 
          next.first_idx, next.current_idx, permute(next.permutation, next.first_idx, next.first_idx - 1)});
      }
      for(int i = next.current_idx; i < weight.size(); i++){
        queue.push({weight[next.permutation[i]] - weight[next.permutation[i - 1]], i - 1, i + 1, permute(next.permutation, i - 1, i)});
      }
    }
    return result;
  }

  std::vector<int> indices(std::vector<int> permutation){
    std::vector<int> indices;
    std::vector<int> result;
    for(int i = 0; i < permutation.size(); i++){
      indices.push_back(i);
    }
    for(int i = 0; i < permutation.size(); i++){
      result.push_back(indices[permutation[i]]);
      for(int j = permutation[i] + 1; j < permutation.size(); j++){
        indices[j]--;
      }
    }
    return result;
  }

  size_t getIndex(std::vector<int> permutation){
    std::vector<int> ids = indices(permutation);
    size_t idx = 0;
    size_t factor = 1;
    size_t product = 1; 
    for(int i = permutation.size() - 2; i >= 0; i--){
      idx += ids[i] * product;
      product *= factor;
    }
    return idx;
  }

  int factorial(int n){
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
  }

  struct PermutationTree{
    struct Node{
      std::vector<Node> children;
      HypernodeWeight avg;
    };
    std::vector<
  }

  struct PermutationDataStructure{
    PartitionedHypergraph& phg;
    std::vector<std::pair<PartitionID, HypernodeID>> id_to_index;
    std::vector<HypernodeID> index_to_id;
    std::vector<std::vector<HypernodeID>> nodes_per_permutation;
    void initialize(){
      std::vector<std::vector<size_t>> prefix_sum;
      prefix_sum.resize(phg.initialNumNodes());
      for(HypernodeID hn : phg.nodes()){
        prefix_sum[hn] = std::vector<size_t>(phg.k(), 0);
        prefix_sum[hn][phg.partID(hn)] = 1;
      }
      parallel_prefix_sum(prefix_sum.begin(), prefix_sum.end(), prefix_sum.begin(), [&](std::vector<size_t> s1, std::vector<size_t> s2){
        std::vector<size_t> res;
        for(int i = 0; i < s1.size(); i++){
          res.push_back(s2[i] + s2[i]);
        }
        return res;
      }, std::vector<size_t>(phg.k(), 0));
      for(HypernodeID hn : phg.nodes()){
        id_to_index[hn] = {phg.partID(hn), prefix_sum[hn][phg.partID(hn)] - 1};
        index_to_id[id_to_index[hn].second] = id_to_index[hn].second;
      }

      

    }
    void insert(HypernodeID hn, std::vector<double> relative_weight){
      nodes_per_permutation[getIndex(indices(getPermutation(relative_weight)))].push_back(hn);
    }
  };*/
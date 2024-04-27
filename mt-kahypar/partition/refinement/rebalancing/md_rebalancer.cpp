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

#include <stdlib.h>

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
    std::vector<std::vector<int>> levels;
    /*std::vector<std::vector<pointer_list>> pointerList;*/
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
      /*pointerList.resize(std::pow(2, imbalanced_dimensions.size()));*/
      levels.resize(std::pow(2, imbalanced_dimensions.size()));
    }

    void initialize(std::vector<std::vector<double>> nn, std::vector<std::pair<double,HypernodeID>> nntw){
      nodes_normalized = nn;
      nodes_normalized_total_weight = nntw;
      for(int i = 0; i < dimension; i++){
        imbalanced_dimensions.push_back(i);
      }
      nodes_prio.resize(std::pow(2, imbalanced_dimensions.size()));
      /*pointerList.resize(std::pow(2, imbalanced_dimensions.size()));*/
      current_indices.resize(std::pow(2, imbalanced_dimensions.size()), 0);
      levels.resize(std::pow(2, imbalanced_dimensions.size()));
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
      /*pointerList[index].resize(nodes_normalized.size() + 1);*/
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
      if(nodes_prio[nodes_prio.size() - 1 - index].size() != 0){
        for(HypernodeID hn = 0; hn < nodes_normalized.size(); hn++){
          nodes_prio[index][hn] = nodes_prio[nodes_prio.size() - 1 - index][nodes_normalized.size() - 1 - hn];
        }
      }
      else{
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
      }

      levels[index].resize(nodes_prio[index].size());
      for(HypernodeID hn = nodes_prio[index].size() - 1; hn-- > 0;){
        HypernodeID node = nodes_prio[index][hn].first;
        int level_number = 0;
        double ignored_percentage = 0.5;
        while(std::ceil(ignored_percentage * nodes_normalized.size()) < nodes_normalized_total_weight[node].second + 1){
          level_number++;
          ignored_percentage = (ignored_percentage + 1.0) / 2.0;
        }
        levels[index][hn] = level_number;
      }



      /*
      std::list<std::pair<int,HypernodeID>> current_last;
      levels[index].resize(nodes_prio[index].size());
      for(HypernodeID hn = nodes_prio[index].size() - 1; hn-- > 0;){
        HypernodeID node = nodes_prio[index][hn].first;
        int level_number = 0;
        double ignored_percentage = 0.5;
        while(std::ceil(ignored_percentage * nodes_normalized.size()) < nodes_normalized_total_weight[node].second + 1){
          level_number++;
          ignored_percentage = (ignored_percentage + 1.0) / 2.0;
        }
        levels[index][hn] = level_number;
        int idx_in_current = 0;
        HypernodeID pointer;
        bool exists_smaller = false;
        for(auto x : current_last){
          if(x.first < level_number){
            exists_smaller =true;
            pointer = x.second;
            continue;
          }
          pointerList[index][hn + 1].next_per_level.push_back(x);
        }
        if(exists_smaller && (pointerList[index][hn + 1].next_per_level.size() == 0 || pointerList[index][hn + 1].next_per_level[0].first != level_number)){
          pointerList[index][hn + 1].next_per_level.insert(pointerList[index][hn + 1].next_per_level.begin(), {level_number, pointer});
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


      for(HypernodeID hn = 0; hn < nodes_prio[index].size() + 1; hn++){
        for(int i = 0; i < pointerList[index][hn].next_per_level.size(); i++){
          std::cout << pointerList[index][hn].next_per_level[i].first << " " << pointerList[index][hn].next_per_level[i].second << "\n";
        }
        std::cout << "end\n";
        
      }*/

      

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
      int max_level = level;
      for(HypernodeID hn = current_indices[index]; hn < nodes_prio[index].size(); hn++){
        max_level = std::max(max_level, levels[index][hn]);
        if(levels[index][hn] <= level){
          current_indices[index] = hn + 1;
          return {nodes_prio[index][hn].first, max_level};
        }
      }
      current_indices[index] = 0;
      return {0, -1};

      /*pointer_list current = pointerList[index][current_indices[index]];    
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
      return {nodes_prio[index][current.next_per_level[start].second - 1].first, lowest_above};*/
    }

    /*struct pointer_list getentry(std::vector<bool> imbalanced, HypernodeID index){
      return pointerList[get_index(imbalanced)][index];
    }*/

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
    /*std::cout << "limits\n";
    for(PartitionID p = 0; p < limits.size(); p++){
      for(int d = 0; d < dimension; d++){
        std::cout << limits[p][d] << "\n";
        if(limits[p][d] < 0.0){
          return {false, packing};
        }
      }
    }*/
    return {true, packing};
  }

  std::pair<bool, std::vector<std::vector<std::pair<PartitionID, HypernodeID>>>> bin_packing_best_fit(std::vector<std::vector<double>> limits, std::vector<std::pair<double, std::pair<PartitionID, HypernodeID>>> nodes, std::vector<std::vector<std::vector<double>>> weights){
    std::vector<std::vector<std::pair<PartitionID, HypernodeID>>> packing(limits.size());
    for(HypernodeID hn = 0; hn < nodes.size(); hn++){
      double min_bg = std::numeric_limits<double>::max();
      PartitionID max_p = -1;
      for(PartitionID p = 0; p < limits.size(); p++){
        std::vector<double> new_limit = minus(limits[p], weights[nodes[hn].second.first][nodes[hn].second.second]);
        bool too_big = false;
        double this_bg = 0.0;
        for(int d = 0; d < dimension; d++){
          if(new_limit[d] <= 0.0){
            too_big = true;
          }
          this_bg += new_limit[d];

        }
        if(too_big){
          continue;
        }
        if(min_bg > this_bg){
          max_p = p;
          min_bg = this_bg;
        }
      }
      if(max_p == -1){
        return {false, packing};
      }
      limits[max_p] = minus(limits[max_p], weights[nodes[hn].second.first][nodes[hn].second.second]);
      packing[max_p].push_back(nodes[hn].second);
    }
    for(PartitionID p = 0; p < limits.size(); p++){
      for(int d = 0; d < dimension; d++){
        if(limits[p][d] < 0.0){
          return {false, packing};
        }
      }
    }
    return {true, packing};
  }





  


  std::pair<bool, std::vector<std::vector<std::pair<PartitionID, HypernodeID>>>> bin_packing_random(std::vector<std::vector<double>> limits, std::vector<std::pair<double, std::pair<PartitionID, HypernodeID>>> nodes, std::vector<std::vector<std::vector<double>>> weights){
    std::vector<std::vector<std::pair<PartitionID, HypernodeID>>> packing(limits.size());
    std::cout << "binpack\n";
    HypernodeID number_random_assertions = std::rand() % nodes.size();
    for(HypernodeID hn = 0; hn < number_random_assertions; hn++){
      PartitionID to = -1;
      for(PartitionID p = 0; p < limits.size(); p++){
        to = std::rand() % limits.size();
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
        break;
      }
      if(to = -1){
        return {false, packing};
      }
      limits[to] = minus(limits[to], weights[nodes[hn].second.first][nodes[hn].second.second]);
      packing[to].push_back(nodes[hn].second);
    }
    for(HypernodeID hn = number_random_assertions; hn < nodes.size(); hn++){
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
    for(PartitionID p = 0; p < limits.size(); p++){
      for(int d = 0; d < dimension; d++){
        if(limits[p][d] < 0.0){
          return {false, packing};
        }
      }
    }
    return {true, packing};
  }


  std::pair<bool, std::vector<std::vector<std::pair<PartitionID, HypernodeID>>>> bin_packing_totally_random(std::vector<std::vector<double>> limits, std::vector<std::pair<double, std::pair<PartitionID, HypernodeID>>> nodes, std::vector<std::vector<std::vector<double>>> weights){
    std::vector<std::vector<std::pair<PartitionID, HypernodeID>>> packing(limits.size());
    std::cout << "binpack\n";
    for(HypernodeID hn = 0; hn < nodes.size(); hn++){
      PartitionID to = -1;
      for(PartitionID p = 0; p < limits.size(); p++){
        to = std::rand() % limits.size();
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
        break;
      }
      if(to = -1){
        return {false, packing};
      }
      limits[to] = minus(limits[to], weights[nodes[hn].second.first][nodes[hn].second.second]);
      packing[to].push_back(nodes[hn].second);
    }
    for(PartitionID p = 0; p < limits.size(); p++){
      for(int d = 0; d < dimension; d++){
        if(limits[p][d] < 0.0){
          return {false, packing};
        }
      }
    }
    return {true, packing};
  }


  /*std::pair<bool, std::vector<std::vector<std::pair<PartitionID, HypernodeID>>>> bin_packing_random(std::vector<std::vector<double>> limits, 
    std::vector<std::pair<double, std::pair<PartitionID, HypernodeID>>> nodes, std::vector<std::vector<std::vector<double>>> weights){
      std::vector<std::vector<std::pair<PartitionID, HypernodeID>>> packing(limits.size());
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
  }*/
  template <typename GraphAndGainTypes>
  bool MDRebalancer<GraphAndGainTypes>::greedyRefiner(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                       vec<vec<Move>>* moves_by_part,
                                                       vec<Move>* moves_linear,
                                                       Metrics& best_metrics, 
                                                       std::vector<HypernodeID> nodes, bool balance){                                                        
                                                   
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph); 
    bool use_max_norm = false;
    //calculate min part weights
    std::vector<HypernodeWeight> min_part_weights(phg.k());
    /*for(PartitionID p = 0; p < phg.k(); p++){
      for(int d = 0; d < dimension; d++){
        min_part_weights[p].weights[d] = 0std::floor(_context.partition.max_part_weights[p].weights[d] * (1.0 - 5.0 * _context.partition.epsilon[0]));
        std::cout << _context.partition.max_part_weights[p].weights[d] << " " << phg.partWeight(p).weights[d] << "\n";
      }
      std::cout << "---------\n";      
    }*/

    std::vector<HypernodeWeight> max_part_weights_modified(phg.k());
    for(PartitionID p = 0; p < phg.k(); p++){
      max_part_weights_modified[p] = 0.9999 * _context.partition.max_part_weights[p];
    }



    auto positive_rebalancing_move = [&](Move_Internal move){
      return move.balance < 0.0 || move.balance <= 0.0 && move.gain > 0;
    };

    auto positive_refinement_move = [&](Move_Internal move){
      return move.gain > 0 || move.gain == 0 && move.balance < 0.0;
    };
    auto is_positive_move = [&](Move_Internal move){
      return balance && positive_rebalancing_move(move) || !balance && positive_refinement_move(move);
    };

    auto weighed_imbalance = [&](){
      double res = 0.0;
      for(PartitionID p = 0; p < phg.k(); p++){
        for(int d = 0; d < dimension; d++){
          res += std::max(0.0, (phg.partWeight(p).weights[d] - _context.partition.max_part_weights[p].weights[d]) * _context.partition.max_part_weights_inv[p][d]);
        }
      }
      return res;
    };

    auto max_imbalance = [&](){
      double res = 0.0;
      for(PartitionID p = 0; p < phg.k(); p++){
        for(int d = 0; d < dimension; d++){
          res = std::max(res, (phg.partWeight(p).weights[d] - _context.partition.max_part_weights[p].weights[d]) * _context.partition.max_part_weights_inv[p][d]);
        }
      }
      return res;
    };
    
    //std::cout << "initial ib: " << weighed_imbalance() << "\n";                                                   
    std::vector<HypernodeID> id_to_index(phg.initialNumNodes(), nodes.size());
    for(HypernodeID h = 0; h < nodes.size(); h++){
      /*if(phg.partID(nodes[h]) == -1){
        std::cout << "error\n\n\n\n\n";
      }*/
      id_to_index[nodes[h]] = h;
    }
    auto objective_delta = [&](const SynchronizedEdgeUpdate& sync_update) {
      _gain.computeDeltaForHyperedge(sync_update);
    };                                                    
    auto horizontal_balance_gain = [&](const PartitionedHypergraph& phg, HypernodeID node, PartitionID from, 
        PartitionID to, const HypernodeWeight min_from = HypernodeWeight(0), const HypernodeWeight min_to = HypernodeWeight(0)){
      double gain = 0.0;
      if(from == to){
        return 0.0;
      }
      for(int i = 0; i < dimension; i++){
        int32_t to_excess = std::max(0, std::min(phg.nodeWeight(node).weights[i], phg.partWeight(to).weights[i] + 
        phg.nodeWeight(node).weights[i] - max_part_weights_modified[to].weights[i]));
        int32_t from_excess = std::max(0, std::min(phg.nodeWeight(node).weights[i], phg.partWeight(from).weights[i] - max_part_weights_modified[from].weights[i]));

        /*int32_t to_undergain = std::max(0, std::min(phg.nodeWeight(node).weights[i], min_to.weights[i] - phg.partWeight(to).weights[i]));
        int32_t from_underloss = std::max(0, std::min(phg.nodeWeight(node).weights[i], min_from.weights[i] - phg.partWeight(from).weights[i] + phg.nodeWeight(node).weights[i]));*/
        if(to_excess == from_excess){
          continue;
        }
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
      /*auto vertical_balance_gain = [&](const PartitionedHypergraph& phg, HypernodeID node, PartitionID from, PartitionID to){
        return vertical_imbalance(phg.partWeight(from) - phg.nodeWeight(node), from) - vertical_imbalance(phg.partWeight(from), from) 
        + vertical_imbalance(phg.partWeight(to) + phg.nodeWeight(node), to) - vertical_imbalance(phg.partWeight(to), to);
      };*/

      auto max_norm_balance_gain = [&](const PartitionedHypergraph& phg, HypernodeID node, PartitionID from, PartitionID to){
        double max = 0.0;
        for(PartitionID p = 0; p < phg.k(); p++){
          for(int d = 0; d < dimension; d++){
            if(phg.partWeight(p).weights[d] / _context.partition.max_part_weights[p].weights[d] > max){
              max = phg.partWeight(p).weights[d] / _context.partition.max_part_weights[p].weights[d];
            }
          }
        }
        double max_after = 0.0;
        for(PartitionID p = 0; p < phg.k(); p++){
          for(int d = 0; d < dimension; d++){
            if(p == from && (phg.partWeight(p).weights[d] - phg.nodeWeight(node).weights[d]) / _context.partition.max_part_weights[p].weights[d] > max_after){
              max_after = (phg.partWeight(p).weights[d] - phg.nodeWeight(node).weights[d]) / _context.partition.max_part_weights[p].weights[d];
            }
            else if(p == to && (phg.partWeight(p).weights[d] + phg.nodeWeight(node).weights[d]) / _context.partition.max_part_weights[p].weights[d] > max_after){
              max_after = (phg.partWeight(p).weights[d] + phg.nodeWeight(node).weights[d]) / _context.partition.max_part_weights[p].weights[d];
            }
            else if(p != from && p != to && phg.partWeight(p).weights[d] / _context.partition.max_part_weights[p].weights[d] > max_after){
              max_after = phg.partWeight(p).weights[d] / _context.partition.max_part_weights[p].weights[d];
            }
          }
        }
        return max_after - max;
      };

      auto is_positive_move_max_norm = [&](const PartitionedHypergraph& phg, HypernodeID node, PartitionID from, PartitionID to, Gain gain){
        return max_norm_balance_gain(phg,node,from,to) < 0.0 || max_norm_balance_gain(phg,node,from,to) == 0.0 && horizontal_balance_gain(phg,node,from,to) < 0.0 
          || max_norm_balance_gain(phg,node,from,to) == 0.0 && horizontal_balance_gain(phg,node,from,to) == 0.0 && gain > 0;
      };

      auto balance_gain = [&](const PartitionedHypergraph& phg, HypernodeID node, PartitionID from, PartitionID to, HypernodeWeight from_min, HypernodeWeight to_min){
        return use_max_norm ? max_norm_balance_gain(phg, node, from, to)  : horizontal_balance_gain(phg, node, from, to);
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
       
      std::vector<tbb::concurrent_vector<std::pair<HypernodeID, HypernodeID>>> nodes_sorted;
      nodes_sorted.resize(mt_kahypar::dimension);
      for(HypernodeID h = 0; h < nodes.size(); h++){
        HypernodeID hn = nodes[h];
        for(int i = 0; i < dimension; i++){
          ASSERT(hn < phg.initialNumNodes());
          nodes_sorted[i].push_back({phg.nodeWeight(hn).weights[i], h});
        }
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
      /*for(int d = 0; d < dimension; d++){
        for(PartitionID p = 0; p < phg.k(); p++){
          indices[i][j] = getBorder(i, exceed[j]);
        }    
      }*/     
      auto get_max_dimension = [&](HypernodeID hn){
        int dim = 0;
        double mw = 0.0;
        for(int d = 0; d < dimension; d++){
          if(_context.partition.max_part_weights_inv[phg.partID(hn)][d] * phg.nodeWeight(hn).weights[d] > mw){
            mw = _context.partition.max_part_weights_inv[phg.partID(hn)][d] * phg.nodeWeight(hn).weights[d];
            dim = d;
          }
        }
        return dim;
      };


      auto get_max_move = [&](HypernodeID node){
        PartitionID p_max = -1;
        Move_Internal max_move = {std::numeric_limits<double>::max(), 0, 0.0};
        for(PartitionID p = 0; p < phg.k(); p++){
          Move_Internal move = {0.0, _gain_cache.gain(node, phg.partID(node), p), balance_gain(phg, node, phg.partID(node), p, min_part_weights[phg.partID(node)], min_part_weights[p])};
          //std::cout << balance_gain(phg, node, phg.partID(node), p, min_part_weights[phg.partID(node)], min_part_weights[p]) << "\n";
          move.recomputeBalance();
          if(is_positive_move(move) && (move.gain_and_balance < max_move.gain_and_balance || p == -1)){
            max_move = move;
            p_max = p;
          }
        }
        return std::pair<PartitionID, Move_Internal>(p_max, max_move);
      };

      Gain quality = metrics::quality(phg, _context);
      std::vector<std::vector<AddressablePQ<HypernodeID, double>>> queues(phg.k());
      std::vector<std::vector<std::vector<HypernodeID>>> index_to_id_partition_dimension(phg.k());
      std::vector<HypernodeID> id_to_index_partition_dimension(phg.initialNumNodes(), phg.initialNumNodes());
      for(PartitionID p = 0; p < phg.k(); p++){
        index_to_id_partition_dimension[p].resize(dimension);
        queues[p].resize(dimension);
      }
      for(HypernodeID hn : nodes){
        id_to_index_partition_dimension[hn] = index_to_id_partition_dimension[phg.partID(hn)][get_max_dimension(hn)].size();
        index_to_id_partition_dimension[phg.partID(hn)][get_max_dimension(hn)].push_back(hn);
      }
      //std::vector<AddressablePQ<PartitionID, Move_internal>> queue_per_node(nodes.size());
      for(PartitionID p = 0; p < phg.k(); p++){
        for(int d = 0; d < dimension; d++){
          queues[p][d].setSize(index_to_id_partition_dimension[p][d].size());
        }
      }
      for(HypernodeID hn : nodes){        
        if(get_max_move(hn).first != -1){
          queues[phg.partID(hn)][get_max_dimension(hn)].insert({id_to_index_partition_dimension[hn], get_max_move(hn).second.gain_and_balance});
        }
      }

      PartitionID imbalanced = 0;
      for(PartitionID p = 0; p < phg.k(); p++){
        if(phg.partWeight(p) > _context.partition.max_part_weights[p]){
          imbalanced++;
        }
      }

      int counter = 0;
      int other_counter = 0;
      const int UPDATE_FREQUENCY = std::ceil(_context.partition.update_frequency * std::sqrt(nodes.size()));
      bool horizontal_balance_used = true;
      std::vector<HypernodeID> min_affected_node;
      for(int d = 0; d < dimension; d++){
        min_affected_node.push_back(nodes.size());
      }
    Gain local_attributed_gain = 0;
    double bag = 0.0;

    /*int i = 1;
      while(i*2 < queue.top_moves.v.size()){
        for(int j = 0; j < i; j++){
          std::cout <<  "(" << queue.top_moves.gains_and_balances[queue.top_moves.v[i - 1 + j]] << "," << queue.top_moves.in_use[queue.top_moves.v[i - 1 + j]] << ") ";
          for(int k = 0; k < queue.queues_per_node[queue.top_moves.v[i - 1 + j]].v.size(); k++){
            std::cout << "(" << queue.queues_per_node[queue.top_moves.v[i - 1 + j]].gains_and_balances[queue.queues_per_node[queue.top_moves.v[i - 1 + j]].v[k]].gain_and_balance << ","
            << queue.queues_per_node[queue.top_moves.v[i - 1 + j]].gains_and_balances[queue.queues_per_node[queue.top_moves.v[i - 1 + j]].v[k]].gain << "," <<
            queue.queues_per_node[queue.top_moves.v[i - 1 + j]].gains_and_balances[queue.queues_per_node[queue.top_moves.v[i - 1 + j]].v[k]].balance << ")";
          }
          std::cout << ") ";
        }
        std::cout << "\n";
        i *= 2;
      }*/
    std::vector<HypernodeWeight> old_part_weights;
    for(PartitionID p = 0; p < phg.k(); p++){
      old_part_weights.push_back(phg.partWeight(p));
    }


    auto queuesEmpty = [&](){
      for(PartitionID p = 0; p < phg.k(); p++){
        for(int d = 0; d < dimension; d++){
          if(!queues[p][d].isEmpty()){
            return false;
          }
        }
      }
      return true;
    };


    std::vector<HypernodeWeight> max_part_weight_hard;
    for(PartitionID p = 0; p < phg.k(); p++){
      max_part_weight_hard.push_back(HypernodeWeight(std::numeric_limits<uint32_t>::max()));
    }

    std::vector<bool> moved(phg.initialNumNodes(), false);
    
    //std::cout << "beforequeue\n";
    bool test_bool = false;
    while(!queuesEmpty() && (balance && (imbalanced != 0) || !balance && max_imbalance() < 0.1)){
      std::pair<PartitionID, Move_Internal> max_move;
      max_move.first = -1;
      std::pair<HypernodeID, double> max_node = {0, std::numeric_limits<double>::max()};
      for(PartitionID p = 0; p < phg.k(); p++){

        for(int d = 0; d < dimension; d++){
          //std::cout << queues[p][d].isEmpty() << "\n";
          if(!queues[p][d].isEmpty()){
            HypernodeID this_node = index_to_id_partition_dimension[p][d][queues[p][d].getMax().first];
            std::pair<PartitionID, Move_Internal> move = get_max_move(this_node);
            double this_prio = move.first == -1 ? queues[p][d].getMax().second :
             std::min(move.second.gain_and_balance, queues[p][d].getMax().second);
            if(this_prio < max_node.second){
              max_node = {this_node, this_prio};
              max_move = move;
            }
            /*std::cout << queues[p][d].getMax().second << " " << max_node.second << " " 
              << get_max_move(index_to_id_partition_dimension[p][d][queues[p][d].getMax().first]).second.gain << " " << 
              get_max_move(index_to_id_partition_dimension[p][d][queues[p][d].getMax().first]).second.balance <<
              " " << get_max_move(index_to_id_partition_dimension[p][d][queues[p][d].getMax().first]).second.gain_and_balance << "\n";*/
          }
          
          /*if(!queues[p][d].isEmpty() && max_node.second > queues[p][d].getMax().second){
            
            max_node = queues[p][d].getMax();
            max_node.first = index_to_id_partition_dimension[p][d][max_node.first];
          }*/
        }
        /*std::cout << "node: " << max_node.first << " " << max_node.second << "\n";*/
      }
      //std::pair<PartitionID, Move_internal> max_move = get_max_move(max_node.first);
      /*ASSERT(max_move.first < nodes.size());
      ASSERT(phg.partID(max_move.first) != -1);
      ASSERT(max_move.second.first != -1);*/
      if(max_move.first == -1){
        /*std::cout << "-1";
        std::cout << queues[phg.partID(max_node.first)][get_max_dimension(max_node.first)].pq.size() << " " << queues[phg.partID(max_node.first)][get_max_dimension(max_node.first)].isEmpty() << "\n";*/
        queues[phg.partID(max_node.first)][get_max_dimension(max_node.first)].deleteMax();
        /*std::cout << queues[phg.partID(max_node.first)][get_max_dimension(max_node.first)].pq.size() << "\n";*/
        other_counter++;
      }
      else if(max_node.second < max_move.second.gain_and_balance){
        if(is_positive_move(max_move.second)){
          queues[phg.partID(max_node.first)][get_max_dimension(max_node.first)].changeKey({id_to_index_partition_dimension[max_node.first], max_move.second.gain_and_balance});
        }
        else{
          queues[phg.partID(max_node.first)][get_max_dimension(max_node.first)].deleteMax();
        }
        
        other_counter++;
      }
      else{
        ASSERT(phg.partID(max_node.first) != -1);
        Move move = {phg.partID(max_node.first), max_move.first, max_node.first, max_move.second.gain};
        bag += max_move.second.gain_and_balance;



        //sanity check
        
        /*for(PartitionID p = 0; p < phg.k(); p++){
          if(!queues[p][0].isEmpty()){
            std::cout << "queue: " << index_to_id_partition_dimension[p][0][queues[p][0].getMax().first] << " " << queues[p][0].getMax().second << "\n";
          }
          if(!queues[p][1].isEmpty()){
            std::cout << "queue: " << index_to_id_partition_dimension[p][1][queues[p][1].getMax().first] << " " << queues[p][1].getMax().second << "\n";
          }
          std::cout << "real: " << max_gain[p] << "\n";
        }*/


        /*HypernodeID hmax = -1;
        PartitionID p = -1;
        std::vector<double> max_gain(phg.k(), std::numeric_limits<double>::max());
        for(PartitionID p = 0; p < phg.k(); p++){
          for(HypernodeID hn : nodes){
            Move_internal m = get_max_move(hn).second;
            m.balance = balance_gain(phg, hn, phg.partID(hn), p, min_part_weights[phg.partID(hn)], min_part_weights[p]);
            m.gain = _gain_cache.gain(hn, phg.partID(hn), p);
            m.recomputeBalance();
            if(m.gain_and_balance < max_gain[phg.partID(hn)] && m.is_positive_move()){
              max_gain[phg.partID(hn)] = m.gain_and_balance;
              hmax = hn;
            }
          }
        }*/






        queues[phg.partID(max_node.first)][get_max_dimension(max_node.first)].deleteMax();
        counter++;      
        imbalanced += (phg.partWeight(move.from) - phg.nodeWeight(move.node) > _context.partition.max_part_weights[move.from])
          - (phg.partWeight(move.from) > _context.partition.max_part_weights[move.from])
          - (phg.partWeight(move.to) > _context.partition.max_part_weights[move.to])
          + (phg.partWeight(move.to) + phg.nodeWeight(move.node) > _context.partition.max_part_weights[move.to]);
        std::vector<HyperedgeID> edges_with_gain_change;
        std::vector<HypernodeID> changed_nodes;
        changed_nodes.push_back(move.node);
        if(moves_linear != NULL){
          moves_linear->push_back(move);
        }
        


        
        
        /*for(PartitionID p = 0; p < phg.k(); p++){
          for(int d = 0; d <dimension; d++){
            std::cout << phg.partWeight(p).weights[d] << " " << _context.partition.max_part_weights[p].weights[d] << "\n";
          }
        }*/
        phg.changeNodePart(_gain_cache, move.node, move.from, move.to, HypernodeWeight(true), []{}, [&](const SynchronizedEdgeUpdate& sync_update) {
                        local_attributed_gain += (sync_update.pin_count_in_to_part_after == 1 ? sync_update.edge_weight : 0) +
                        (sync_update.pin_count_in_from_part_after == 0 ? -sync_update.edge_weight : 0);
                        if (!PartitionedHypergraph::is_graph && GainCache::triggersDeltaGainUpdate(sync_update)) {
                          edges_with_gain_change.push_back(sync_update.he);
                        }
                      });
        moved[move.node] = true;
        for(HyperedgeID& he : edges_with_gain_change){
          for(HypernodeID hn : phg.pins(he)){
            changed_nodes.push_back(hn);
          }
        }
        if constexpr (PartitionedHypergraph::is_graph) {
          for (const auto e : phg.incidentEdges(move.node)) {
            HypernodeID v = phg.edgeTarget(e);
            changed_nodes.push_back(v);
          }
        }


        changed_nodes.push_back(move.node);
        if(counter % UPDATE_FREQUENCY == 0/*queuesEmpty()*/){
          
          


          for(int i = 0; i < dimension; i++){          
            HypernodeID min_affected_node = nodes_sorted[i].size() - 1;
            for(PartitionID p = 0; p < phg.k(); p++){
              if(phg.partWeight(p).weights[i] < _context.partition.max_part_weights[p].weights[i] && phg.partWeight(p).weights[i] < old_part_weights[p].weights[i]){
                int prior_diff = _context.partition.max_part_weights[p].weights[i] - old_part_weights[p].weights[i];
                while(min_affected_node >= 0 && nodes_sorted[i][min_affected_node].first > prior_diff){
                  changed_nodes.push_back(nodes[nodes_sorted[i][min_affected_node].second]);
                  if(min_affected_node == 0) break;
                  min_affected_node--;
                }
              }
              else if(phg.partWeight(p).weights[i] > _context.partition.max_part_weights[p].weights[i] && phg.partWeight(p).weights[i] > old_part_weights[p].weights[i]){
                int prior_diff = old_part_weights[p].weights[i] -  _context.partition.max_part_weights[p].weights[i];
                while(min_affected_node >= 0 && nodes_sorted[i][min_affected_node].first > prior_diff){
                  changed_nodes.push_back(nodes[nodes_sorted[i][min_affected_node].second]);
                  if(min_affected_node == 0) break;
                  min_affected_node--;
                }
              }
            }            
          } 
          for(PartitionID p = 0; p < phg.k(); p++){
            old_part_weights[p] = phg.partWeight(p);
          }       
        }
        for(size_t i = 0; i < changed_nodes.size(); i++){
          if(id_to_index_partition_dimension[changed_nodes[i]] != phg.initialNumNodes() && !moved[changed_nodes[i]]){
            std::pair<PartitionID, Move_Internal> best_move = get_max_move(changed_nodes[i]);
            if(!is_positive_move(best_move.second)){
              queues[phg.partID(changed_nodes[i])][get_max_dimension(changed_nodes[i])].invalidate(id_to_index_partition_dimension[changed_nodes[i]]);
            }
            else{
              queues[phg.partID(changed_nodes[i])][get_max_dimension(changed_nodes[i])].changeKey({id_to_index_partition_dimension[changed_nodes[i]], best_move.second.gain_and_balance});
            }
          }          
        }
        bool wrong = false;
        /*for(HypernodeID hn : phg.nodes()){
          if(get_max_move(hn).first != -1 && !moved[hn]){
            if(get_max_move(hn).second.gain_and_balance < queues[phg.partID(hn)][get_max_dimension(hn)].items[id_to_index_partition_dimension[hn]] && !wrong){
              wrong = true;
              std::cout << counter << " " << get_max_move(hn).first << " " << get_max_move(hn).second.gain_and_balance << " " << queues[phg.partID(hn)][get_max_dimension(hn)].items[id_to_index_partition_dimension[hn]] << " " << move.from << " " << move.to << " " << phg.partID(hn) << " firsterror " << get_max_move(hn).second.balance << " " 
              << balance_gain(phg, hn, phg.partID(hn), get_max_move(hn).first, min_part_weights[phg.partID(hn)], min_part_weights[get_max_move(hn).first]) << "\n";
              for(int d = 0; d < dimension; d++){
                std::cout << phg.nodeWeight(hn).weights[d] << "\n";
              }
            }
          }
        }*/
      }
      //if(test_bool){std::cout << "test2\n";}
      if(balance && queuesEmpty() && !use_max_norm && imbalanced != 0){
        std::cout << "new fallback\n";
        //use_max_norm=true;
        test_bool = true;
        double total_overload = 0.0;

        /*for(PartitionID p = 0; p < phg.k(); p++){
          for(HypernodeID hn : phg.nodes()){
            if(!moved[hn] && balance_gain(phg, hn, phg.partID(hn), p, min_part_weights[phg.partID(hn)], min_part_weights[p]) < 0.0){
              //std::cout << "error\n";
            }
          }
        }*/


        for(PartitionID p = 0; p < phg.k(); p++){
          for(int d = 0; d < dimension; d++){
            total_overload += std::max(0.0, static_cast<double>(phg.partWeight(p).weights[d] - _context.partition.max_part_weights[p].weights[d]) / static_cast<double>(_context.partition.max_part_weights[p].weights[d]));
            std::cout << _context.partition.max_part_weights[p].weights[d] << " " << phg.partWeight(p).weights[d] << "\n";
          }
          std::cout << "---------\n";  
        }
        /*for(PartitionID p = 0; p < phg.k(); p++){
          for(int d = 0; d < dimension; d++){
            max_part_weights_modified[p].weights[d] = std::floor((1.0 - _context.partition.epsilon[0] / 8.0) * max_part_weights_modified[p].weights[d]);
            //std::cout << max_part_weights_modified[p].weights[d] << "\n"; 
          }
        }*/

        for(HypernodeID hn : nodes){
          if(moved[hn]) continue;
          std::pair<PartitionID, Move_Internal> best_move = get_max_move(hn);
          if(best_move.first == -1){            
            queues[phg.partID(hn)][get_max_dimension(hn)].invalidate(id_to_index_partition_dimension[hn]);
          }
          else{
            queues[phg.partID(hn)][get_max_dimension(hn)].changeKey({id_to_index_partition_dimension[hn], best_move.second.gain_and_balance});
          }          
        }
      }                                                  
    }
    


    /*for(HypernodeID h = 0; h < nodes.size(); h++){
      for(PartitionID p = 0; p < phg.k(); p++){
        if(queue.queues_per_node[h].gains_and_balances[p].gain != _gain_cache.gain(nodes[h], phg.partID(nodes[h]), p)){
          std::cout << "error!!!!!!!!\n\n\n\n\n";
        }
      }
    }*/

    /*int i = 1;
      while(i*2 < queue.top_moves.v.size()){
        for(int j = 0; j < i; j++){
          std::cout <<  "(" << queue.top_moves.gains_and_balances[queue.top_moves.v[i - 1 + j]] << "," << queue.top_moves.in_use[queue.top_moves.v[i - 1 + j]] << ") ";
          for(int k = 0; k < queue.queues_per_node[queue.top_moves.v[i - 1 + j]].v.size(); k++){
            std::cout << "(" << queue.queues_per_node[queue.top_moves.v[i - 1 + j]].gains_and_balances[queue.queues_per_node[queue.top_moves.v[i - 1 + j]].v[k]].gain_and_balance << ","
            << queue.queues_per_node[queue.top_moves.v[i - 1 + j]].gains_and_balances[queue.queues_per_node[queue.top_moves.v[i - 1 + j]].v[k]].gain << "," <<
            queue.queues_per_node[queue.top_moves.v[i - 1 + j]].gains_and_balances[queue.queues_per_node[queue.top_moves.v[i - 1 + j]].v[k]].balance << ")";
          }
          std::cout << ") ";
        }
        std::cout << "\n";
        i *= 2;
      }*/


    std::cout << counter << " " << other_counter << " " << local_attributed_gain << " " << imbalanced << " " << queuesEmpty() << "\n";
    return imbalanced == 0;
  }

  struct S_Calculator{
    double factor_lowerBound = 0.0;
    double factor_maxbound;
    double factor_upperBound;
    double factor = 1.0;
    double search_threshold = 0.001;
    PartitionID k;
    std::vector<std::vector<double>> max_remaining_weight;
    std::vector<int> current_level_per_partition;
    std::vector<std::vector<std::vector<double>>> normalized_weights;
    //the second part of the pair is its rank in sorted total weights
    std::vector<std::vector<std::pair<double,HypernodeID>>> totalWeight;
    std::vector<std::vector<double>> normalized_partition_weights;
    std::vector<std::vector<double>> normalized_perfect_weights;
    std::vector<std::vector<double>> virtual_weight;

    std::vector<struct hypernodes_ordered> nodes_ordered;
    bool last = false;



    void initialize(std::vector<std::vector<std::vector<double>>> &nws, 
        std::vector<std::vector<std::pair<double,HypernodeID>>> &tW,
        std::vector<std::vector<double>> &npws,
        std::vector<std::vector<double>> &nperfws, double epsilon){
      factor_maxbound = 1.0 + epsilon;
      factor_upperBound = factor_maxbound;
      normalized_weights = nws;
      totalWeight = tW;
      normalized_partition_weights = npws;
      normalized_perfect_weights = nperfws;
      k = normalized_partition_weights.size();
      nodes_ordered.resize(k);
      virtual_weight.resize(k);
      max_remaining_weight.resize(k);
      current_level_per_partition.resize(k,0);
      for(PartitionID p = 0; p < k; p++){
        nodes_ordered[p].initialize(normalized_weights[p], totalWeight[p]);
      }
      current_level_per_partition.resize(k, 0);
    }

    bool hasNext(){
      return !last && (factor_upperBound - factor_lowerBound > search_threshold || factor_lowerBound == 0.0 && factor != 0.0);
    }

    void adjustFactor(bool success){
      if(factor == 0.0){
        last = true;
      }
      if(success){
        factor_lowerBound = factor;
      }
      else{
        factor_upperBound = factor;
      }
      if(factor_lowerBound == 0.0){
        factor = 2.0 * factor_upperBound - factor_maxbound;
      }
      else{
        factor = (factor_upperBound + factor_lowerBound) / 2.0;
      }
      if(factor_lowerBound == 0.0 && factor < search_threshold){
        factor = 0.0;
      }
    }
    
    std::vector<std::vector<bool>> calculate_S(bool success){ 
      

      auto smaller_equal = [](std::vector<double> a, std::vector<double> b){
        for(size_t i = 0; i < a.size(); i++){
          if(a[i] > b[i] + 0.00000000001){
            return false;
          }
        }
        return true;
      };
      std::vector<std::vector<double>> max_remaining_weight(k);
      for(PartitionID p = 0; p < k; p++){
        max_remaining_weight[p].resize(dimension);
        for(int d = 0; d < dimension; d++){
          max_remaining_weight[p][d] = normalized_perfect_weights[p][d] * factor;
        }
      }
      std::vector<std::vector<bool>> extracted(k);
      for(PartitionID p = 0; p < k; p++){
        extracted[p].resize(normalized_weights[p].size(), false);
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
          bool possible_to_balance = true;
          while(!smaller_equal(this_run_weight, max_remaining_weight[p])){
            double max = std::numeric_limits<double>::min();
            double min = std::numeric_limits<double>::max();
            for(int d = 0; d < dimension; d++){
              imbalanced[d] = this_run_weight[d] > max_remaining_weight[p][d];
              max = std::max(max, this_run_weight[d]);
              min = std::min(min, this_run_weight[d]);
            }
            for(int d = 0; d < dimension; d++){
              imbalanced[d] = this_run_weight[d] > (max + min) / 2.0;
            }
            std::pair<HypernodeID,int> next;
            do{
              next = nodes_ordered[p].getNext(imbalanced, level + 1);
            }while(next.second != -1 && used[next.first]);
            if(next.second == -1){
              possible_to_balance = false;
              break;
            }
            used[next.first] = true;
            highest_would_use = std::max(highest_would_use, next.second);
            for(int d = 0; d < dimension; d++){
              this_run_weight[d] -= normalized_weights[p][next.first][d];
            }
          }
              
          nodes_ordered[p].reset();
          if(!possible_to_balance){
            level++;
            continue;
          }
          if(smaller_equal(this_run_weight, max_remaining_weight[p])){
            double quality = 0.0;
            for(int d = 0; d < dimension; d++){
              quality += std::pow(max_remaining_weight[p][d] - this_run_weight[d], 2);
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
              
              
          level = std::max(level + 1, highest_would_use);
        }
        current_level_per_partition[p] = level;
      }
      return extracted;
    }
  

    std::vector<std::vector<double>> get_virtual_weight(){
      return virtual_weight;
    }
  };
  

  template <typename GraphAndGainTypes>
  bool MDRebalancer<GraphAndGainTypes>::refineInternal(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                       vec<vec<Move>>* moves_by_part,
                                                       vec<Move>* moves_linear,
                                                       Metrics& best_metrics) {                                       
    vec<Move> moves;
    auto start = std::chrono::high_resolution_clock::now();                                                
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
    resizeDataStructuresForCurrentK();
    // This function is passed as lambda to the changeNodePart function and used
    // to calculate the "real" delta of a move (in terms of the used objective function).
    auto objective_delta = [&](const SynchronizedEdgeUpdate& sync_update) {
      _gain.computeDeltaForHyperedge(sync_update);
    };


    std::vector<HypernodeID> all_nodes;
      for(HypernodeID hn : phg.nodes()){
        if(phg.partID(hn) != -1){
          all_nodes.push_back(hn);
        }
      }
    for(int i = 0; i < 2; i++){
      vec<HypernodeWeight> initial_part_weights;
      for(PartitionID p = 0; p < phg.k(); p++){
        initial_part_weights.push_back(phg.partWeight(p));
      }
      vec<Move> refine_moves_before;
      Gain before_quality = metrics::quality(phg, _context);
      std::cout << "beforelp " << before_quality << "\n";
      std::cout << _context.partition.allowed_imbalance_refine << "\n";
      labelPropagation(hypergraph, moves_by_part, &refine_moves_before, best_metrics, all_nodes, _context.partition.allowed_imbalance_refine, false);
      //greedyRefiner(hypergraph, moves_by_part, &moves, best_metrics, all_nodes, false);
      vec<Move> rebalance_moves_before;
      std::cout << "before greedy " << metrics::quality(phg, _context) << "\n";
      greedyRefiner(hypergraph, moves_by_part, &rebalance_moves_before, best_metrics, all_nodes, true);
      std::cout << "after greedy " << metrics::quality(phg, _context) << "\n";

      std::vector<std::vector<std::pair<MoveID,Move>>> refine_by_part(phg.k()); 
      std::vector<HypernodeID> refine_idx(phg.k(), 0);
      vec<vec<Move>> rebalance_by_part(phg.k());
      std::vector<HypernodeID> rebalance_idx(phg.k(), 0);
      std::vector<Move> final_moves;
      std::vector<PartitionID> moved_from_refine(phg.initialNumNodes(), -1);
      std::vector<int> last_move_of_node(phg.initialNumNodes());
      std::vector<int> first_move_of_node(phg.initialNumNodes(), -1);
      std::vector<bool> doubly_moved(phg.initialNumNodes(), false);
      vec<Move> refine_moves;
      vec<Move> rebalance_moves;
      for(int idx = 0; idx < refine_moves_before.size(); idx++){
        Move move = refine_moves_before[idx];
        last_move_of_node[move.node] = idx;
        if(first_move_of_node[move.node] == -1){
          first_move_of_node[move.node] = idx;
        }
      }
      for(int idx = 0; idx < refine_moves_before.size(); idx++){
        Move move = refine_moves_before[idx];
        if(last_move_of_node[move.node] == idx){
          ASSERT(first_move_of_node[move.node] != -1);
          PartitionID from = refine_moves_before[first_move_of_node[move.node]].from;
          refine_moves.push_back({from, move.to, move.node, 0});
        }
      }
      for(int idx = 0; idx < rebalance_moves_before.size(); idx++){
        ASSERT(phg.partID(rebalance_moves_before[idx].node) == rebalance_moves_before[idx].to);
        ASSERT(rebalance_moves_before[idx].from != -1);
        ASSERT(rebalance_moves_before[idx].to != -1);
        rebalance_by_part[rebalance_moves_before[idx].from].push_back(rebalance_moves_before[idx]);
      }
      
      /*for(int idx = 0; idx < rebalance_moves_before.size(); idx++){
        Move move = rebalance_moves_before[idx];
        if(last_move_of_node[move.node] >= refine_moves_before.size()){
          std::cout << "errorhere\n";
        }
        last_move_of_node[move.node] = idx + refine_moves_before.size();
        if(moved_from_refine[idx] != -1){
          doubly_moved[idx] = true;
          rebalance_moves_before[idx].from = moved_from_refine[idx];
          if(rebalance_moves_before[idx].from != rebalance_moves_before[idx].to){
            rebalance_moves.push_back(rebalance_moves_before[idx]);
          }
        }
        else{
          rebalance_moves.push_back(rebalance_moves_before[idx]);
        }
      }
      for(int idx = 0; idx < refine_moves_before.size(); idx++){
        Move move = refine_moves_before[idx];
        if(last_move_of_node[move.node] == idx){
          move.from = moved_from_refine[move.node];
          refine_moves.push_back(move);
        }
      }
      for(int i = 0; i < refine_moves.size() + rebalance_moves.size(); i++){
        for(int j = i + 1; j < rebalance_moves.size() + refine_moves.size(); j++){
          HypernodeID n1 = i < refine_moves.size() ? refine_moves[i].node : rebalance_moves[i - refine_moves.size()].node;
          HypernodeID n2 = j < refine_moves.size() ? refine_moves[j].node : rebalance_moves[j - refine_moves.size()].node;
          ASSERT(n1 != n2);
        }
      }*/
      /*for(int idx = 0; idx < rebalance_moves_before.size(); idx++){
        Move move = rebalance_moves_before[rebalance_moves_before.size() - 1 - idx];
        phg.changeNodePart(_gain_cache, move.node, move.to, move.from, HypernodeWeight(true), []{}, [&](const SynchronizedEdgeUpdate& sync_update) {});
      }
      for(int idx = 0; idx < refine_moves.size(); idx++){
        Move move = refine_moves[refine_moves.size() - 1 - idx];
        phg.changeNodePart(_gain_cache, move.node, move.to, move.from, HypernodeWeight(true), []{}, [&](const SynchronizedEdgeUpdate& sync_update) {});
      }*/
      //std::cout << "rolled back moves " << metrics::quality(phg, _context) << "\n";
      vec<Move> move_order;
      interleaveMoveSequenceWithRebalancingMoves(hypergraph, initial_part_weights,
                                                            _context.partition.max_part_weights,
                                                            refine_moves,
                                                            rebalance_by_part,
                                                            move_order);
      //std::cout << "movessize" << move_order.size() << "\n";
      for(int idx = 0; idx < move_order.size(); idx++){
        Move move = move_order[idx];
        ASSERT(move.from != -1);
        ASSERT(move.to != -1);
        ASSERT(phg.partID(move.node) == move.to);
        phg.changeNodePart(_gain_cache, move.node, move.to, move.from, HypernodeWeight(true), []{}, [&](const SynchronizedEdgeUpdate& sync_update) {});
      }
      
      //std::cout << "rolled back " << metrics::quality(phg, _context) << "\n";
      Gain best_quality = 0;
      Gain current_quality = 0;
      int bq_idx = 0;
      for(int idx = 0; idx < move_order.size(); idx++){
        Move move = move_order[idx];
        current_quality += _gain_cache.gain(move.node, move.from,move.to);
        phg.changeNodePart(_gain_cache, move.node, move.from, move.to, HypernodeWeight(true), []{}, [&](const SynchronizedEdgeUpdate& sync_update) {});
        /*if(metrics::isBalanced(phg, _context)){
          std::cout << "index: " << idx << "\n";
          std::cout << "quality: " << current_quality << "\n";
        }*/
        if(metrics::isBalanced(phg, _context) && current_quality > best_quality){
          best_quality = current_quality;
          bq_idx = idx + 1;
        }
      }
      for(int idx = bq_idx; idx < move_order.size(); idx++){
        Move move = move_order[idx];
        phg.changeNodePart(_gain_cache, move.node, move.to, move.from, HypernodeWeight(true), []{}, [&](const SynchronizedEdgeUpdate& sync_update) {});
      }
      std::cout << "after rollback " << metrics::quality(phg, _context) << "\n";
      //std::cout << "bestidx " << bq_idx << "\n";
      ASSERT(before_quality >= metrics::quality(phg, _context));

      /*for(int i = 0; i < refine_moves.size(); i++){
        Move move = refine_moves[i];
        refine_by_part[move.from].push_back({i, move});
        if(move.node == 65047){std::cout << "rm\n";}
        //std::cout << moved_from_refine[move.node] << "\n";
        phg.changeNodePart(_gain_cache, move.node, move.to, move.from, HypernodeWeight(true), []{}, [&](const SynchronizedEdgeUpdate& sync_update) {});
      }
      for(Move move : rebalance_moves){
        rebalance_by_part[move.from].push_back(move);
        if(move.node == 65047){std::cout << "rb\n";}
        phg.changeNodePart(_gain_cache, move.node, move.to, move.from, HypernodeWeight(true), []{}, [&](const SynchronizedEdgeUpdate& sync_update) {});
      }
      Gain quality_change = 0;
      Gain best_qc = 0;
      int max_idx = 0;
      int not_at_end = 0;
      if(refine_moves.size() > 0){not_at_end++;}
      for(PartitionID p = 0; p < phg.k(); p++){
        if(rebalance_by_part[p].size() > 0){
          not_at_end++;
        }
      }
      while(true){
        PartitionID move_p = -1;
        Move move;
        int index = 0;
        int loop = 0;
        for(PartitionID p = 0; p < phg.k(); p++){
          if(!(phg.partWeight(p) <= _context.partition.max_part_weights[p])){
            move_p = p;
            break;
          }
        }
        PartitionID next_refine_move = -1;
        MoveID position = std::numeric_limits<MoveID>::max();
        for(PartitionID p = 0; p < phg.k(); p++){
          if(refine_idx[p] < refine_by_part[p].size() && position > refine_by_part[p][refine_idx[p]]){
            position = refine_by_part[p][refine_idx[p]];
            next_refine_move = p;
          }
        }
        if(move_p != -1){
          if(refine_idx[move_p] < refine_by_part[move_p].size()){
            move = refine_by_part[move_p][refine_idx[move_p]].second;
            refine_idx[move_p]++;
          }
          else{
            ASSERT(rebalance_idx[move_p] < rebalance_by_part[move_p].size());
            move = rebalance_by_part[move_p][rebalance_idx[move_p]];
            
            loop = 1;
            rebalance_idx[move_p]++;
            index = rebalance_idx[move_p];
            if(rebalance_idx[move_p] >= rebalance_by_part[move_p].size()){
              not_at_end--;
            }
          }
          
        }
        else if(refine_idx < refine_moves.size()){
          move = refine_moves[refine_idx];
          index = refine_idx;
          loop = 2;
          refine_idx++;
          if(refine_idx >= refine_moves.size()){
            not_at_end--;
          }
        }
        else{
          loop = 3;
          PartitionID not_empty = -1;
          for(PartitionID p = 0; p < phg.k(); p++){
            if(rebalance_by_part[p].size() > rebalance_idx[p]){
              not_empty = p;
              break;
            }
          }
          if(not_empty == -1){break;}
          move = rebalance_by_part[not_empty][rebalance_idx[not_empty]];
          index = rebalance_idx[not_empty];
          
          rebalance_idx[not_empty]++;
        }
        quality_change += _gain_cache.gain(move.node, move.from, move.to);
        if(phg.partID(move.node) != move.from || move.node == 65047){
          std::cout << move.node << " " << move.from << " " << phg.partID(move.node) << " " << index << " " << loop << " " << move_p << "\n";
        }
        ASSERT(phg.partID(move.node) == move.from, move.from);
        phg.changeNodePart(_gain_cache, move.node, move.from, move.to, HypernodeWeight(true), []{}, [&](const SynchronizedEdgeUpdate& sync_update) {});
        final_moves.push_back(move);
        if(quality_change > best_qc && metrics::isBalanced(phg, _context)){
          best_qc = quality_change;
          max_idx = final_moves.size();
        }
      }
      while(true){}
      for(int m = max_idx; m < final_moves.size(); m++){
        Move move = final_moves[m];
        phg.changeNodePart(_gain_cache, move.node, move.to, move.from, 
          HypernodeWeight(true), []{}, [&](const SynchronizedEdgeUpdate& sync_update) {});
      }*/

    }
    labelPropagation(hypergraph, moves_by_part, &moves, best_metrics, all_nodes, 0.0, false);
    
    // If partition is imbalanced, rebalancer is activated
    if ( !metrics::isBalanced(phg, _context) ) {
      DBG << "Starting multi-dimensional rebalancer";  // only printed if debug=true in header
      _gain.reset();
      if(!greedyRefiner(hypergraph, moves_by_part, &moves, best_metrics, all_nodes, true)){
        for(PartitionID p = 0; p < phg.k(); p++){
          for(int d = 0; d < dimension; d++){
            std::cout << phg.partWeight(p).weights[d] << " " << _context.partition.max_part_weights[p].weights[d] << "\n";
          }
        }
        std::vector<HypernodeID> L;
        const double L_threshold = _context.partition.node_threshold;

        //weights of the partition just taking V/L into account
        std::vector<HypernodeWeight> partition_weight_included;

        std::vector<std::vector<std::vector<double>>> normalized_weights(phg.k());
        //total Weight of each node (organized per partition)
        std::vector<std::vector<std::pair<double,HypernodeID>>> totalWeight(phg.k());
        std::vector<std::vector<std::pair<double,HypernodeID>>> totalWeight_ordered(phg.k());
        std::vector<std::vector<double>> normalized_partition_weights(phg.k());
        std::vector<std::vector<double>> normalized_perfect_weights(phg.k());
        std::vector<std::pair<PartitionID, bool>> imbalances_per_partition(phg.k());
        std::vector<struct hypernodes_ordered> nodes_ordered(phg.k());
        std::vector<HypernodeID> map_id_to_index(phg.initialNumNodes());
        std::vector<std::vector<HypernodeID>> map_index_to_id(phg.k());

        std::vector<HypernodeWeight> part_weight_without_L(phg.k());
      

      
        for(int p = 0; p < phg.k(); p++){
          part_weight_without_L[p] = HypernodeWeight(0);
          normalized_partition_weights[p].resize(dimension);
          normalized_perfect_weights[p].resize(dimension);
          for(int d = 0; d < dimension; d++){
            normalized_partition_weights[p][d] = 0.0;
            partition_weight_included.push_back(HypernodeWeight(0));
            normalized_perfect_weights[p][d] = _context.partition.perfect_balance_part_weights[p].weights[d] * _context.partition.max_part_weights_inv[p][d];
          }
        }

        for(HypernodeID hn : phg.nodes()){
          if(phg.partID(hn) == -1){
            continue;
          }
          PartitionID part = phg.partID(hn);
          double total_weight = 0.0;
          std::vector<double> weight;
          double max_weight = 0.0;
          for(int d = 0; d < dimension; d++){
            double w = phg.nodeWeight(hn).weights[d] * _context.partition.max_part_weights_inv[part][d];
            total_weight += w;
            weight.push_back(w);
            max_weight = std::max(max_weight, w);         
          }
          if(max_weight > L_threshold){
            part_weight_without_L[part] = part_weight_without_L[part] + phg.nodeWeight(hn);
            map_id_to_index[hn] = normalized_weights[part].size();
            map_index_to_id[part].push_back(hn);
            normalized_weights[part].push_back(weight);
            totalWeight_ordered[part].push_back({total_weight, map_id_to_index[hn]});
            partition_weight_included[part] += phg.nodeWeight(hn);
          }
          else{
            L.push_back(hn);
          }
          
        }

        bool overloaded = false;
        for(PartitionID p = 0; p < phg.k(); p++){
          overloaded = overloaded || partition_weight_included[p] > _context.partition.max_part_weights[p];
        }

        if(overloaded){

          for(PartitionID p = 0; p < phg.k(); p++){
            for(int d = 0 ; d < dimension; d++){
              normalized_partition_weights[p][d] = part_weight_without_L[p].weights[d] * _context.partition.max_part_weights_inv[p][d];
            }
          }

          for(PartitionID p = 0; p < phg.k(); p++){
            std::sort(totalWeight_ordered[p].begin(), totalWeight_ordered[p].end());
            totalWeight[p].resize(totalWeight_ordered[p].size());
            for(HypernodeID hn = 0; hn < totalWeight_ordered[p].size(); hn++){
              totalWeight[p][(totalWeight_ordered[p][hn]).second] = {(totalWeight_ordered[p][hn]).first, hn};
            }
          }

          struct S_Calculator s_calculator;
          s_calculator.initialize(normalized_weights, totalWeight, normalized_partition_weights, 
            normalized_perfect_weights, _context.partition.epsilon[0]);
          bool success = true;
          std::vector<std::vector<std::pair<PartitionID, HypernodeID>>> current_result(phg.k());
          while(s_calculator.hasNext()){
            std::vector<std::vector<bool>> extracted = s_calculator.calculate_S(success);
            std::vector<std::pair<double, std::pair<PartitionID,HypernodeID>>> bp_nodes;
            std::vector<std::pair<double, std::pair<PartitionID,HypernodeID>>> bp_nodes_absolute;
            std::vector<std::vector<double>> limits(phg.k());
            for(PartitionID p = 0; p < phg.k(); p++){
              for(int d = 0; d < dimension; d++){
                limits[p].push_back(1.0);
              }
              for(HypernodeID hn = 0; hn < extracted[p].size(); hn++){
                if(extracted[p][hn]){
                  double min = std::numeric_limits<double>::max();
                  double max = std::numeric_limits<double>::min();
                  for(int d = 0; d < dimension; d++){
                    min = std::min(min, normalized_weights[p][hn][d]);
                    max = std::max(max, normalized_weights[p][hn][d]);
                    
                  }
                  bp_nodes.push_back({(max - min) / max, {p, hn}});
                  bp_nodes_absolute.push_back({(max - min), {p, hn}});
                }
                else{
                  for(int d = 0; d < dimension; d++){
                    limits[p][d] -= normalized_weights[p][hn][d];
                  }
                }
              }
            }
            std::sort(bp_nodes.begin(), bp_nodes.end(), [&](auto a, auto b){
              return a.first > b.first;
            });
            std::pair<bool, std::vector<std::vector<std::pair<PartitionID, HypernodeID>>>> result = bin_packing(limits, bp_nodes, normalized_weights);
            if(result.first){
              current_result = result.second;
            }
            else{
              std::sort(bp_nodes.begin(), bp_nodes.end(), [&](auto a, auto b){
                return totalWeight[a.second.first][a.second.second].first > totalWeight[b.second.first][b.second.second].first;
              });
              result = bin_packing_best_fit(limits, bp_nodes, normalized_weights);
              if(result.first){
              current_result = result.second;
            }
            }
            /*else{
              for(int i = 0; i < 20; i++){
                result = bin_packing_random(limits, bp_nodes, normalized_weights);
                if(result.first){
                  current_result = result.second;
                  break;
                }
              }
              for(int i = 0; i < 20; i++){
                result = bin_packing_totally_random(limits, bp_nodes, normalized_weights);
                if(result.first){
                  current_result = result.second;
                  break;
                }
              }
            }*/
            success = result.first;
            //std::cout << s_calculator.factor << " " << s_calculator.factor_lowerBound << " " << s_calculator.factor_upperBound << " " << success << "\n";
            s_calculator.adjustFactor(success);
              
          }
          Gain local_attributed_gain = 0;
          for(PartitionID p = 0; p < phg.k(); p++){
            for(HypernodeID hn = 0; hn < current_result[p].size(); hn++){
              HypernodeID node = map_index_to_id[current_result[p][hn].first][current_result[p][hn].second];
              /*phg.changeNodePart(node, phg.partID(node), p, objective_delta);*/
              if(phg.partID(node) == p){
                continue;
              }
              moves.push_back({phg.partID(node), p, node, _gain_cache.gain(node, phg.partID(node), p)});
              
              phg.changeNodePart(_gain_cache, node, phg.partID(node), p, HypernodeWeight(true), []{}, [&](const SynchronizedEdgeUpdate& sync_update) {
                          local_attributed_gain += (sync_update.pin_count_in_to_part_after == 1 ? sync_update.edge_weight : 0) +
                          (sync_update.pin_count_in_from_part_after == 0 ? -sync_update.edge_weight : 0);
                        });
            }
          }
          //std::cout << "fallback result: " << local_attributed_gain << "\n";
        }
        //return greedyRefiner(hypergraph, moves_by_part, moves_linear, best_metrics, L);
        

        



        

        
      
      }
      labelPropagation(hypergraph, moves_by_part, &moves, best_metrics, all_nodes, 0.0, false);
      // TODO: rebalancing logic goes here
    
    //bool res = greedyRefiner(hypergraph, moves_by_part, &moves, best_metrics, all_nodes);

    if (moves_by_part != nullptr) {
      moves_by_part->resize(_context.partition.k);
      for (auto& direction : *moves_by_part) direction.clear();
      for (size_t i = 0; i < moves.size(); ++i) {
        (*moves_by_part)[(moves[i]).from].push_back(moves[i]);
      }
    } /*else if (moves_linear != nullptr) {
      moves_linear->clear();
      moves_linear->reserve(num_moves_performed);
      for (size_t i = 0; i < num_moves_performed; ++i) {
        moves_linear->push_back(_moves[i]);
      }
    }*/



      /*if (moves_by_part != nullptr) {
        moves_by_part->resize(_context.partition.k);
        for (auto& direction : *moves_by_part) direction.clear();
        // TODO: ignore for now, implementation necessary to support unconstrained refinement (compare advanced_rebalancer)
      } else if (moves_linear != nullptr) {
        moves_linear->clear();
      // TODO: ignore for now, implementation necessary to support unconstrained refinement (compare advanced_rebalancer)
      }*/

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
    
      
    }
    return true;
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




  template <typename GraphAndGainTypes>
  bool MDRebalancer<GraphAndGainTypes>::labelPropagation(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                       vec<vec<Move>>* moves_by_part,
                                                       vec<Move>* moves_linear,
                                                       Metrics& best_metrics, 
                                                       std::vector<HypernodeID> nodes, double allowed_imbalance, bool balance){

    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
    std::vector<HypernodeWeight> max_part_weights(phg.k());
    auto get_imbalances = [&](){
      std::vector<double> ib(phg.k(), 0.0);
      for(PartitionID p = 0; p < phg.k(); p++){
        for(int i = 0; i < dimension; i++){
          ib[i] = std::max(ib[i], (phg.partWeight(p).weights[i] - _context.partition.max_part_weights[p].weights[i]) *
           _context.partition.max_part_weights_inv[p][i]);
        }
      }
      return ib;
    };      
    

    std::vector<double> imbalances = get_imbalances();
    for(PartitionID p = 0; p < phg.k(); p++){
      for(int i = 0; i < dimension; i++){
        max_part_weights[p].weights[i] = std::ceil((1.0 + allowed_imbalance + imbalances[i]) * _context.partition.max_part_weights[p].weights[i]); 
      }
    }

    if(balance){
      for(PartitionID p = 0; p < phg.k(); p++){
        for(int i = 0; i < dimension; i++){
          max_part_weights[p].weights[i] = _context.partition.max_part_weights[p].weights[i]; 
        }
      }
    }
    else{
      for(PartitionID p = 0; p < phg.k(); p++){
        for(int i = 0; i < dimension; i++){
          max_part_weights[p].weights[i] = std::max(max_part_weights[p].weights[i], _context.partition.max_part_weights[p].weights[i]); 
        }
      }
    }

    auto is_smaller = [&](int n, int scalar, HypernodeWeight x, HypernodeWeight y, HypernodeWeight z){
      for(int i = 0; i < n; i++){
        if(scalar * x.weights[i] + y.weights[i] > z.weights[i]){
          return false;
        }
      }
      return true;
    };


    auto BetterBalanceKWay = [&](HypernodeWeight vwgt, 
        int a1, PartitionID p1, 
        int a2, PartitionID p2)
    { 
      double nrm1=0.0, nrm2=0.0, max1=0.0, max2=0.0;

      for (int i = 0; i < dimension; i++) {
        double tmp = (phg.partWeight(p1).weights[i] + a1 * vwgt.weights[i] - max_part_weights[p1].weights[i]) / max_part_weights[p1].weights[i];
        //printf("BB: %d %+.4f ", (int)i, (float)tmp);
        nrm1 += tmp*tmp;
        max1 = (tmp > max1 ? tmp : max1);

        tmp = (phg.partWeight(p2).weights[i] + a2 * vwgt.weights[i] - max_part_weights[p2].weights[i]) / max_part_weights[p2].weights[i];
        //printf("%+.4f ", (float)tmp);
        nrm2 += tmp*tmp;
        max2 = (tmp > max2 ? tmp : max2);

      }

      if (max2 < max1)
        return 1;

      if (max2 == max1 && nrm2 < nrm1)
        return 1;

      return 0;
    };




    AddressablePQ<HypernodeID,double> pq;
    std::vector<std::pair<HypernodeID,HypernodeID>> id_ed(phg.initialNumNodes());
    std::vector<std::vector<HypernodeID>> num_neighbors(phg.initialNumNodes());
    std::vector<PartitionID> nnp(phg.initialNumNodes(), 0);
    for(HypernodeID hn : phg.nodes()){
      id_ed[hn] = std::pair<HypernodeID,HypernodeID>(0,0);
      num_neighbors[hn].resize(phg.k(), 0);
      for(HyperedgeID he : phg.incidentEdges(hn)){
        for(HypernodeID h : phg.pins(he)){
          if(h == hn) continue;
          if(phg.partID(h) == phg.partID(hn)){
            id_ed[hn].first += phg.edgeWeight(he);
          }
          else{
            id_ed[hn].second += phg.edgeWeight(he);
            if(num_neighbors[hn][phg.partID(h)]== 0){
              nnp[hn]++;
            }
          }
          num_neighbors[hn][phg.partID(h)] += phg.edgeWeight(he);
        }
      }
    }

    //std::cout << num_neighbors[14942][0] << " " << id_ed[14942].first << " thiss\n\n";
    
    std::vector<HypernodeID> idx_in_bnd(phg.initialNumNodes(), 0);
    int niter = 15;
    for(int run = 0; run < niter; run++){
      std::vector<HypernodeID> bnd;
      for(HypernodeID hn : phg.nodes()){
        if(balance){
          if(id_ed[hn].second > 0){
            bnd.push_back(hn);
          }
        }
        else{
          if(id_ed[hn].second >= id_ed[hn].first){
            bnd.push_back(hn);
          }
        }
      }
      pq.setSize(phg.initialNumNodes());
      for(HypernodeID hn : bnd){
        double gain = (nnp[hn] > 0 ? (id_ed[hn].second / sqrt(nnp[hn])) : 0.0) - id_ed[hn].first;
        pq.insert({hn, -gain});
      }
      Gain quality_change = 0;
      int i = 0;
      while(!pq.isEmpty()){
        i++;
        HypernodeID hn = pq.deleteMax().first;
        /*if (!balance) {          
          if (id_ed[hn].first > 0 && 
              
              !ivecaxpygez(dimension, -1, vwgt+i*ncon, pwgts+from*ncon, minpwgts+from*ncon))
            continue;   
        }
        else { /* OMODE_BALANCE */
          /*if (!ivecaxpygez(ncon, -1, vwgt+i*ncon, pwgts+from*ncon, minpwgts+from*ncon)) 
            continue;   
        }*/
        PartitionID from = phg.partID(hn);
        PartitionID max_to = -1;
        Gain max_gain = 0; 
        if (!balance) {
                 
          for(PartitionID k = phg.k() - 1; k >=0; k--){
            Gain gain = num_neighbors[hn][k] - num_neighbors[hn][phg.partID(hn)];
            if(max_to != -1 && (num_neighbors[hn][k] > num_neighbors[hn][max_to] && is_smaller(dimension, 1, phg.nodeWeight(hn), phg.partWeight(k), max_part_weights[max_to])
              || num_neighbors[hn][k] == num_neighbors[hn][max_to] && BetterBalanceKWay(phg.nodeWeight(hn), 1, max_to, 1, k))){
              max_to = k;
              max_gain = gain;
            }
            
            if (max_to == -1 && gain >= 0 && is_smaller(dimension, 1, phg.nodeWeight(hn), phg.partWeight(k), max_part_weights[k])){
              max_to = k;
              max_gain = gain;
            }
        
          }
          if (max_to == -1)
            continue;  
          if (!(max_gain > 0 || (max_gain == 0 && BetterBalanceKWay(phg.nodeWeight(hn), -1, phg.partID(hn), +1, max_to)))){
            continue;
          }
            
            
        }
        else {  /* OMODE_BALANCE */
          /*for (k=myrinfo->nnbrs-1; k>=0; k--) {
            if (!safetos[to=mynbrs[k].pid])
              continue;
            if (ivecaxpylez(ncon, 1, vwgt+i*ncon, pwgts+to*ncon, maxpwgts+to*ncon) || 
                BetterBalanceKWay(ncon, vwgt+i*ncon, ubfactors,
                    -1, pwgts+from*ncon, pijbm+from*ncon,
                    +1, pwgts+to*ncon, pijbm+to*ncon))
              break;
          }
          if (k < 0)
            continue;  /* break out if you did not find a candidate */

          /*cto = to;
          for (j=k-1; j>=0; j--) {
            if (!safetos[to=mynbrs[j].pid])
              continue;
            if (BetterBalanceKWay(ncon, vwgt+i*ncon, ubfactors, 
                    1, pwgts+cto*ncon, pijbm+cto*ncon,
                    1, pwgts+to*ncon, pijbm+to*ncon)) {
              k   = j;
              cto = to;
            }
          }
          to = cto;

          if (mynbrs[k].ed-myrinfo->id < 0 &&
              !BetterBalanceKWay(ncon, vwgt+i*ncon, ubfactors,
                    -1, pwgts+from*ncon, pijbm+from*ncon,
                    +1, pwgts+to*ncon, pijbm+to*ncon))
            continue;*/
        }



        /*=====================================================================
        * If we got here, we can now move the vertex from 'from' to 'to' 
        *======================================================================*/
        Move move = {phg.partID(hn), max_to, hn, max_gain};
        moves_linear->push_back(move);
        phg.changeNodePart(_gain_cache, hn, phg.partID(hn), max_to);
        quality_change += num_neighbors[hn][from] - num_neighbors[hn][max_to];
        //std::cout << "pc: " << quality_change << " " << hn << " ";


        for(HyperedgeID he : phg.incidentEdges(hn)){
          for(HypernodeID h : phg.pins(he)){
            if(h != hn){
              if(num_neighbors[h][max_to] == 0 && phg.partID(h) != max_to){
                nnp[h]++;
              }
              num_neighbors[h][from] -= phg.edgeWeight(he);
              num_neighbors[h][max_to] += phg.edgeWeight(he);
              if(num_neighbors[h][from] == 0 && phg.partID(h) != from){
                nnp[h]--;
              }
              if(phg.partID(h) == from){
                id_ed[h].first -= phg.edgeWeight(he);
                id_ed[h].second += phg.edgeWeight(he);
              }
              if(phg.partID(h) == max_to){
                id_ed[h].second -= phg.edgeWeight(he);
                id_ed[h].first += phg.edgeWeight(he);
              }
              if(balance){
                if(id_ed[hn].second > 0){
                  pq.changeKey({h, -((nnp[h] > 0 ? (id_ed[h].second / sqrt(nnp[h])) : 0.0) - id_ed[h].first)});
                }
                else{
                  pq.invalidate(h);
                }
              }
              else{
                if(id_ed[hn].second > 0 && id_ed[hn].second >= id_ed[hn].first){
                  pq.changeKey({h, -((nnp[h] > 0 ? (id_ed[h].second / sqrt(nnp[h])) : 0.0) - id_ed[h].first)});
                }
                else{
                  pq.invalidate(h);
                }
              }
            }
          }
        }

        /*IFSET(ctrl->dbglvl, METIS_DBG_MOVEINFO, 
            printf("\t\tMoving %6"PRIDX" to %3"PRIDX". Gain: %4"PRIDX". Cut: %6"PRIDX"\n", 
                i, to, mynbrs[k].ed-myrinfo->id, graph->mincut));

        /* Update ID/ED and BND related information for the moved vertex */
        /*iaxpy(ncon,  1, vwgt+i*ncon, 1, pwgts+to*ncon,   1);
        iaxpy(ncon, -1, vwgt+i*ncon, 1, pwgts+from*ncon, 1);
        UpdateMovedVertexInfoAndBND(i, from, k, to, myrinfo, mynbrs, where, 
            nbnd, bndptr, bndind, bndtype);
        
        /* Update the degrees of adjacent vertices */
        /*for (j=xadj[i]; j<xadj[i+1]; j++) {
          ii = adjncy[j];
          me = where[ii];
          myrinfo = graph->ckrinfo+ii;

          oldnnbrs = myrinfo->nnbrs;

          UpdateAdjacentVertexInfoAndBND(ctrl, ii, xadj[ii+1]-xadj[ii], me, 
              from, to, myrinfo, adjwgt[j], nbnd, bndptr, bndind, bndtype);

          UpdateQueueInfo(queue, vstatus, ii, me, from, to, myrinfo, oldnnbrs, 
              nupd, updptr, updind, bndtype);

          ASSERT(myrinfo->nnbrs <= xadj[ii+1]-xadj[ii]);
        }*/
      }
      if(quality_change==0 && !balance) break;
    }
    

    /*graph->nbnd = nbnd;

    /* Reset the vstatus and associated data structures */
    /*for (i=0; i<nupd; i++) {
      ASSERT(updptr[updind[i]] != -1);
      ASSERT(vstatus[updind[i]] != VPQSTATUS_NOTPRESENT);
      vstatus[updind[i]] = VPQSTATUS_NOTPRESENT;
      updptr[updind[i]]  = -1;
    }

    if (ctrl->dbglvl&METIS_DBG_REFINE) {
       printf("\t[%6"PRIDX" %6"PRIDX"], Bal: %5.3"PRREAL", Nb: %6"PRIDX"."
              " Nmoves: %5"PRIDX", Cut: %6"PRIDX", Vol: %6"PRIDX,
              imin(nparts*ncon, pwgts,1), imax(nparts*ncon, pwgts,1), 
              ComputeLoadImbalance(graph, nparts, pijbm), 
              graph->nbnd, nmoved, graph->mincut, ComputeVolume(graph, where));
       if (ctrl->minconn) 
         printf(", Doms: [%3"PRIDX" %4"PRIDX"]", imax(nparts, nads,1), isum(nparts, nads,1));
       printf("\n");
    }

    if (nmoved == 0 || (omode == OMODE_REFINE && graph->mincut == oldcut))
      break;
  }

  rpqDestroy(queue);
    }  */ 
    return true;     
  }

  template<typename GraphAndGainTypes>
  void MDRebalancer<GraphAndGainTypes>::interleaveMoveSequenceWithRebalancingMoves(
                                                            mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                            const vec<HypernodeWeight>& initialPartWeights,
                                                            const std::vector<HypernodeWeight>& max_part_weights,
                                                            vec<Move>& refinement_moves,
                                                            vec<vec<Move>>& rebalancing_moves_by_part,
                                                            vec<Move>& move_order) {
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
    ASSERT(rebalancing_moves_by_part.size() == static_cast<size_t>(_context.partition.k));
    /*for(PartitionID p = 0; p < phg.k(); p++){
      for(Move move : rebalancing_moves_by_part[p]){
        ASSERT(move.from != -1);
        ASSERT(move.to != -1);
      }
    }*/
    /*HEAVY_REFINEMENT_ASSERT([&] {
      std::set<HypernodeID> moved_nodes;
      for (PartitionID part = 0; part < context.partition.k; ++part) {
        for (const Move& m: rebalancing_moves_by_part[part]) {
          if (m.from != part || m.to != phg.partID(m.node) || moved_nodes.count(m.node) != 0) {
            return false;
          }
          moved_nodes.insert(m.node);
        }
      }
      return true;
    }());*/

    std::vector<int> wasNodeMovedInThisRound(phg.initialNumNodes(), -1);
    for(int idx = 0; idx < refinement_moves.size(); idx++){
      Move move = refinement_moves[idx];
      ASSERT(wasNodeMovedInThisRound[move.node] == -1);
      wasNodeMovedInThisRound[move.node] = idx;
    }
    // Check the rebalancing moves for nodes that are moved twice. Double moves violate the precondition of the global
    // rollback, which requires that each node is moved at most once. Thus we "merge" the moves of any node
    // that is moved twice (e.g., 0 -> 2 -> 1 becomes 0 -> 1)
    for (PartitionID part = 0; part < _context.partition.k; ++part) {
      vec<Move>& moves = rebalancing_moves_by_part[part];
      tbb::parallel_for(UL(0), moves.size(), [&](const size_t i) {
        Move& r_move = moves[i];
        if (r_move.isValid() && wasNodeMovedInThisRound[r_move.node] != -1 && refinement_moves[wasNodeMovedInThisRound[r_move.node]].isValid()) {
          ASSERT(r_move.to == phg.partID(r_move.node));
          Move& first_move = refinement_moves[wasNodeMovedInThisRound[r_move.node]];
          ASSERT(r_move.node == first_move.node && r_move.from == first_move.to);
          if (first_move.from == r_move.to) {
            // if rebalancing undid the move, we simply delete it
            first_move.invalidate();
            r_move.invalidate();
          } else {
            // "merge" the moves
            r_move.from = first_move.from;
            first_move.invalidate();
          }
        }
      }, tbb::static_partitioner());
      /*for(Move m : moves){
        ASSERT(!m.isValid() || m.to == phg.partID(m.node));
      }*/
    }
    /*for(Move m : refinement_moves){
        ASSERT(!m.isValid() || m.to == phg.partID(m.node));
      }*/

    // NOTE: We re-insert invalid rebalancing moves to ensure the gain cache is updated correctly by the global rollback
    // For now we use a sequential implementation, which is probably fast enough (since this is a single scan trough
    // the move sequence). We might replace it with a parallel implementation later.
    vec<HypernodeWeight> current_part_weights = initialPartWeights;
    vec<MoveID> current_rebalancing_move_index(_context.partition.k, 0);
    MoveID next_move_index = 0;

    auto insert_moves_to_balance_part = [&](const PartitionID part) {
      if (current_part_weights[part] > max_part_weights[part]) {
        insertMovesToBalanceBlock(hypergraph, part, max_part_weights, rebalancing_moves_by_part,
                                  next_move_index, current_part_weights, current_rebalancing_move_index, move_order);
      }
    };

    // it might be possible that the initial weights are already imbalanced
    for (PartitionID part = 0; part < _context.partition.k; ++part) {
      insert_moves_to_balance_part(part);
    }

    const MoveID num_moves = refinement_moves.size();
    for (MoveID move_id = 0; move_id < num_moves; ++move_id) {
      const Move& m = refinement_moves[move_id];
      if (m.isValid()) {
        const HypernodeWeight hn_weight = phg.nodeWeight(m.node);
        current_part_weights[m.from] -= hn_weight;
        current_part_weights[m.to] += hn_weight;
        ASSERT(m.from != -1);
        ASSERT(m.to != -1);
        move_order.push_back(m);
        ++next_move_index;
        // insert rebalancing moves if necessary
        insert_moves_to_balance_part(m.to);
      } else {
        // setting moveOfNode to zero is necessary because, after replacing the move sequence,
        // wasNodeMovedInThisRound() could falsely return true otherwise
        
      }
    }

    // append any remaining rebalancing moves (rollback will decide whether to keep them)
    for (PartitionID part = 0; part < _context.partition.k; ++part) {
      while (current_rebalancing_move_index[part] < rebalancing_moves_by_part[part].size()) {
        const MoveID move_index_for_part = current_rebalancing_move_index[part];
        const Move& m = rebalancing_moves_by_part[part][move_index_for_part];
        ++current_rebalancing_move_index[part];
        if(m.from != -1){
          move_order.push_back(m);
        }
        ++next_move_index;
      }
    }
  }

  template<typename GraphAndGainTypes>
  void MDRebalancer<GraphAndGainTypes>::insertMovesToBalanceBlock(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                                        const PartitionID block,
                                                                        const std::vector<HypernodeWeight>& max_part_weights,
                                                                        const vec<vec<Move>>& rebalancing_moves_by_part,
                                                                        MoveID& next_move_index,
                                                                        vec<HypernodeWeight>& current_part_weights,
                                                                        vec<MoveID>& current_rebalancing_move_index,
                                                                        vec<Move>& move_order) {
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
    while (current_part_weights[block] > max_part_weights[block]
            && current_rebalancing_move_index[block] < rebalancing_moves_by_part[block].size()) {
      const MoveID move_index_for_block = current_rebalancing_move_index[block];
      const Move& m = rebalancing_moves_by_part[block][move_index_for_block];
      ++current_rebalancing_move_index[block];
      /*ASSERT(block < _context.partition.k);
      ASSERT(m.from != -1, block);
      ASSERT(m.to != -1);*/
      
      ++next_move_index;
      if (m.isValid()) {
        move_order.push_back(m);
        const HypernodeWeight hn_weight = phg.nodeWeight(m.node);
        current_part_weights[m.from] -= hn_weight;
        current_part_weights[m.to] += hn_weight;

        if (current_part_weights[m.to] > max_part_weights[m.to]) {
          // edge case: it is possible that the rebalancing move itself causes new imbalance -> call recursively
          insertMovesToBalanceBlock(hypergraph, m.to, max_part_weights, rebalancing_moves_by_part,
                                    next_move_index, current_part_weights, current_rebalancing_move_index, move_order);
        }
      }
    }
  }

  /*int ivecaxpygez(int n, idx_t a, idx_t *x, idx_t *y, idx_t *z)
  {
    for (n--; n>=0; n--) {
      if (a*x[n]+y[n] < z[n]) 
        return 0;
    }

    return  1;
  }*/
}




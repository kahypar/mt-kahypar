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
        if(true/*levels[index][hn] <= level*/){
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
    std::cout << "limits\n";
    for(PartitionID p = 0; p < limits.size(); p++){
      for(int d = 0; d < dimension; d++){
        std::cout << limits[p][d] << "\n";
        if(limits[p][d] < 0.0){
          return {false, packing};
        }
      }
    }
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
                                                       std::vector<HypernodeID> nodes){                                                        

    std::cout << "start\n";                                                      
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph); 
    auto weighed_imbalance = [&](){
      double res = 0.0;
      for(PartitionID p = 0; p < phg.k(); p++){
        for(int d = 0; d < dimension; d++){
          res += std::max(0.0, (phg.partWeight(p).weights[d] - _context.partition.max_part_weights[p].weights[d]) * _context.partition.max_part_weights_inv[p][d]);
        }
      }
      return res;
    };
    std::cout << "initial ib: " << weighed_imbalance() << "\n";                                                   
    std::vector<HypernodeID> id_to_index(phg.initialNumNodes(), nodes.size());
    for(HypernodeID h = 0; h < nodes.size(); h++){
      if(phg.partID(nodes[h]) == -1){
        std::cout << "error\n\n\n\n\n";
      }
      id_to_index[nodes[h]] = h;
    }
    auto objective_delta = [&](const SynchronizedEdgeUpdate& sync_update) {
      _gain.computeDeltaForHyperedge(sync_update);
    };                                                    
    auto horizontal_balance_gain = [&](const PartitionedHypergraph& phg, HypernodeID node, PartitionID from, PartitionID to){
      double gain = 0.0;
      if(from == to){
        return 0.0;
      }
      for(int i = 0; i < dimension; i++){
        int32_t to_excess = std::max(0, std::min(phg.nodeWeight(node).weights[i], phg.partWeight(to).weights[i] + 
        phg.nodeWeight(node).weights[i] - _context.partition.max_part_weights[to].weights[i]));
        int32_t from_excess = std::max(0, std::min(phg.nodeWeight(node).weights[i], phg.partWeight(from).weights[i] - _context.partition.max_part_weights[from].weights[i]));
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
      Gain quality = metrics::quality(phg, _context);
      MoveQueue queue;
      queue.initialize(nodes.size(), phg.k());
      for(HypernodeID h = 0; h < nodes.size(); h++){
        HypernodeID hn = nodes[h];
        for(PartitionID p = 0; p < phg.k(); p++){
          if(p == phg.partID(hn)){
            continue;
          }
          Move_internal move = {0.0, _gain_cache.gain(hn, phg.partID(hn), p), balance_gain(phg, hn, phg.partID(hn), p, true)};
          move.recomputeBalance();
          queue.insert_without_updating({h, {p, move}});
        }
      }
      queue.check();
      PartitionID imbalanced = 0;
      for(PartitionID p = 0; p < phg.k(); p++){
        if(phg.partWeight(p) > _context.partition.max_part_weights[p]){
          imbalanced++;
        }
      }

      int counter = 0;
      const int UPDATE_FREQUENCY = 50;
      bool horizontal_balance_used = true;
      std::vector<HypernodeID> min_affected_node;
      for(int d = 0; d < dimension; d++){
        min_affected_node.push_back(nodes.size());
      }
    Gain local_attributed_gain = 0;
    double bag = 0.0;
    while(imbalanced != 0 && !queue.isEmpty()){
      std::pair<HypernodeID, std::pair<PartitionID, Move_internal>> max_move = queue.getMax();
      if(_gain_cache.gain(max_move.first, phg.partID(nodes[max_move.first]), max_move.second.first) != max_move.second.second.gain){
        std::cout << _gain_cache.gain(max_move.first, phg.partID(nodes[max_move.first]), max_move.second.first) << "  " <<max_move.second.second.gain << "\n";
      }
      ASSERT(max_move.first < nodes.size());
      ASSERT(phg.partID(max_move.first) != -1);
      ASSERT(max_move.second.first != -1);
      HypernodeID node = max_move.first;
      PartitionID from = phg.partID(nodes[node]);
      while(from == -1){
        std::cout << "test\n";
      }
      PartitionID to = max_move.second.first;
      double balance = balance_gain(phg, nodes[node], from, to, horizontal_balance_used);
      if(balance > max_move.second.second.balance){
        queue.changeBalance({node, {to, balance}});
        queue.update(node);
      }
      else{
        bag += max_move.second.second.balance;
        queue.deleteMax();
        counter++;
        Move_internal gains = max_move.second.second;
        HypernodeWeight from_old = phg.partWeight(from);
        HypernodeWeight to_old = phg.partWeight(to);
        HypernodeWeight exceed_from_old = exceed[from];
        HypernodeWeight exceed_to_old = exceed[to];
        Move move = {phg.partID(nodes[node]), to, nodes[node], 0};
        imbalanced += (phg.partWeight(from) - phg.nodeWeight(nodes[node]) > _context.partition.max_part_weights[from])
          - (phg.partWeight(from) > _context.partition.max_part_weights[from])
          - (phg.partWeight(to) > _context.partition.max_part_weights[to])
          + (phg.partWeight(to) + phg.nodeWeight(nodes[node]) > _context.partition.max_part_weights[to]);
        while(to == -1){
          std::cout << "error\n";
        }
        std::vector<HyperedgeID> edges_with_gain_change;
        std::vector<HypernodeID> changed_nodes;
        changed_nodes.push_back(nodes[node]);
        phg.changeNodePart(_gain_cache, nodes[node], from, to, HypernodeWeight(true), []{}, [&](const SynchronizedEdgeUpdate& sync_update) {
                        local_attributed_gain += (sync_update.pin_count_in_to_part_after == 1 ? sync_update.edge_weight : 0) +
                        (sync_update.pin_count_in_from_part_after == 0 ? -sync_update.edge_weight : 0);
                        if (!PartitionedHypergraph::is_graph && GainCache::triggersDeltaGainUpdate(sync_update)) {
                          edges_with_gain_change.push_back(sync_update.he);
                        }
                      });
        for(HyperedgeID& he : edges_with_gain_change){
          for(HypernodeID hn : phg.pins(he)){
            changed_nodes.push_back(hn);
          }
        }
        if constexpr (PartitionedHypergraph::is_graph) {
          for (const auto e : phg.incidentEdges(node)) {
            HypernodeID v = phg.edgeTarget(e);
            changed_nodes.push_back(v);
          }
        }
        for(HypernodeID hn : changed_nodes){
          for(PartitionID p = 0; p < phg.k(); p++){
            if(id_to_index[hn] != nodes.size()){
              queue.changeGain({id_to_index[hn], {p, _gain_cache.gain(hn, phg.partID(hn), p)}});
            }
          }
        }
        calcExceed(from);
        calcExceed(to);
        /*tbb::concurrent_vector<HypernodeID> changed_gain_nodes = _gain.getChangedMoves(phg, {move.from, move.to, move.node, 0}, &queue, id_to_index, nodes.size());
        std::vector<HypernodeID> changed_nodes;
        for(HypernodeID h = 0; h < changed_gain_nodes.size(); h++){
          changed_nodes.push_back(id_to_index[changed_gain_nodes[h]]);
        }*/
        /*for(PartitionID p = 0; p < phg.k(); p++){
          queue.addToGain({node, {p, -gains.gain}});
          queue.changeBalance({node, {p, balance_gain(phg, nodes[node], to, p, horizontal_balance_used)}});
        }*/
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
          if(counter % UPDATE_FREQUENCY == 0){
            counter++;
            
            std::cout <<"gain: " << local_attributed_gain << " " << weighed_imbalance() << "\n";
            for(HypernodeID hn = min_affected_node[i]; hn < nodes_sorted[i].size(); hn++){
              for(PartitionID p = 0; p < phg.k(); p++){
                queue.changeBalance({nodes_sorted[i][hn].second, {p, balance_gain(phg, nodes[nodes_sorted[i][hn].second], 
                phg.partID(nodes[nodes_sorted[i][hn].second]), p, horizontal_balance_used)}});
              }
              changed_nodes.push_back(nodes[nodes_sorted[i][hn].second]);
            }               
          }
        }
        for(size_t i = 0; i < changed_nodes.size(); i++){
          if(id_to_index[changed_nodes[i]] != nodes.size()){
            queue.update(id_to_index[changed_nodes[i]]);
          }          
        }
        queue.checkSizes(phg.k());
      }                                                  
    }
    std::cout << counter << " " << local_attributed_gain << "\n";
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

      std::vector<HypernodeID> all_nodes;
      for(HypernodeID hn : phg.nodes()){
        if(phg.partID(hn != -1)){
          all_nodes.push_back(hn);
        }
      }

      if(true/*!greedyRefiner(hypergraph, moves_by_part, moves_linear, best_metrics, all_nodes)*/){
        std::cout << "fb activated";
        std::vector<HypernodeID> L;
        const double L_threshold = _context.partition.node_threshold;

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
            normalized_perfect_weights[p][d] = _context.partition.perfect_balance_part_weights[p].weights[d] * _context.partition.max_part_weights_inv[p][d];
          }
        }

        //determine size of L
        /*std::vector<std::pair<double, HypernodeID>> nodes_max_weight;
        for(HypernodeID hn : phg.nodes()){
          double max_node_weight = 0.0;
          for(int d = 0; d < dimension; d++){
            max_node_weight = std::max(max_node_weight, phg.nodeWeight(hn).weights[d] * _context.partition.max_part_weights_inv[phg.partID(hn)][d]);
          }
          nodes_max_weight.push_back({max_node_weight, hn});
        }
        std::sort(nodes_max_weight.begin(), nodes_max_weight.end());
        for(HypernodeID i = 0; i < nodes_max_weight.size(); i++){

        }*/

        //initialize normalized_weights, totalWeight and normalized_partition_weights
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
            /*for(int d = 0; d < dimension; d++){
              normalized_partition_weights[part][d] += weight[d];
            }*/
            map_id_to_index[hn] = normalized_weights[part].size();
            map_index_to_id[part].push_back(hn);
            normalized_weights[part].push_back(weight);
            totalWeight_ordered[part].push_back({total_weight, map_id_to_index[hn]});
          }
          else{
            L.push_back(hn);
          }
          
        }
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
          for(PartitionID p = 0; p < phg.k(); p++){
            for(int d = 0; d < dimension; d++){
              std::cout << limits[p][d] << "\n";
            }
          }
          std::pair<bool, std::vector<std::vector<std::pair<PartitionID, HypernodeID>>>> result = bin_packing(limits, bp_nodes, normalized_weights);
          std::cout << "success: " << result.first << "\n";
          if(result.first){
            current_result = result.second;
          }
          else{
            std::sort(bp_nodes.begin(), bp_nodes.end(), [&](auto a, auto b){
              return totalWeight[a.second.first][a.second.second].first > totalWeight[b.second.first][b.second.second].first;
            });
            result = bin_packing_best_fit(limits, bp_nodes, normalized_weights);
            if(result.first){
              std::cout << "succeeee\n";
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
          std::cout << s_calculator.factor << " " << s_calculator.factor_lowerBound << " " << s_calculator.factor_upperBound << " " << success << "\n";
          s_calculator.adjustFactor(success);
            
        }
        Gain local_attributed_gain = 0;
        for(PartitionID p = 0; p < phg.k(); p++){
          for(HypernodeID hn = 0; hn < current_result[p].size(); hn++){
            HypernodeID node = map_index_to_id[current_result[p][hn].first][current_result[p][hn].second];
            /*phg.changeNodePart(node, phg.partID(node), p, objective_delta);*/
            phg.changeNodePart(_gain_cache, node, phg.partID(node), p, HypernodeWeight(true), []{}, [&](const SynchronizedEdgeUpdate& sync_update) {
                        local_attributed_gain += (sync_update.pin_count_in_to_part_after == 1 ? sync_update.edge_weight : 0) +
                        (sync_update.pin_count_in_from_part_after == 0 ? -sync_update.edge_weight : 0);
                      });
          }
        }
        std::cout << "fallback result: " << local_attributed_gain << "\n";
        return greedyRefiner(hypergraph, moves_by_part, moves_linear, best_metrics, L);


        



        

        
      
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


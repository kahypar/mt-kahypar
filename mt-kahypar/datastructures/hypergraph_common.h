/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#pragma once

#include <cstdint>
#include <limits>

#include "include/libmtkahypartypes.h"

#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/hardware_topology.h"
#include "mt-kahypar/parallel/tbb_initializer.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/datastructures/array.h"

namespace mt_kahypar {
static constexpr size_t dimension = 2;
struct NodeWeight {
  std::array<int32_t, dimension> weights;

  NodeWeight(const NodeWeight &nw) {
    for(int i = 0; i < dimension; i++){
        weights[i] = nw.weights[i];
    }
  }

  explicit NodeWeight(){
    for(int i = 0; i < dimension; i++){
        weights[i] = 0;
    }
  }

  explicit NodeWeight(bool maxValue){
    if(maxValue){
      for(int i = 0; i < dimension; i++){
        weights[i] = std::numeric_limits<int32_t>::max();
      }
    }
    else{
      for(int i = 0; i < dimension; i++){
        weights[i] = std::numeric_limits<int32_t>::min();
      }
    }
    
  }

  explicit NodeWeight(int value){
    for(int i = 0; i < dimension; i++){
      weights[i] = value;
    }
  }

  explicit NodeWeight(std::array<double, dimension> d){
    for(int i = 0; i < dimension; i++){
      weights[i] = std::floor(d[i]);
    }
  }
  explicit NodeWeight(std::array<double, dimension> d, bool floor){
    for(int i = 0; i < dimension; i++){
      if(floor){
        weights[i] = std::floor(d[i]);
      }
      else{
        weights[i] = std::ceil(d[i]);
      }
    }
  }

  NodeWeight scale(double d) const{
    NodeWeight res;
    for(int i = 0; i < dimension; i++){
      res.weights[i] = std::ceil(static_cast<double>(weights[i]) * d);
    }
    return res;
  }
  /*constexpr NodeWeight(const uint32_t value, bool b){
    for(int i = 0; i < dimension; i++){
      weights[i] = value;
    }
  }*/

  int32_t max(){
    int32_t max = weights[0];
    for(int i = 1; i < dimension; i++){
      max = std::max(max, weights[i]);
    }
    return max;
  }

  int32_t min(){
    int32_t minim = weights[0];
    for(int i = 1; i < dimension; i++){
      minim = std::min(minim, weights[i]);
    }
    return minim;
  }


  NodeWeight operator =(const NodeWeight &nw){
    for(int i = 0; i < dimension; i++){
      weights[i] = nw.weights[i];
    }
    return *this;
  }

  bool operator >(const NodeWeight &nw) const{
    for(int i = 0; i < dimension; i++){
      if(weights[i] > nw.weights[i]){
        return true;
      }
    }
    return false;
  }

  bool operator <(const NodeWeight &nw) const{
    for(int i = 0; i < dimension; i++){
      if(weights[i] >= nw.weights[i]){
        return false;
      }
    }
    return true;
  }

  
  NodeWeight cutToZero() const{
    NodeWeight res;
    for(int i = 0; i < dimension; i++){
      res.weights[i] = std::max(0, weights[i]);
    }
    return res;
  }

  NodeWeight(uint32_t weight){
    for(int i = 0; i < dimension; i++){
      weights[i] = 0;
    }
  }

  NodeWeight operator +(const NodeWeight &ew) const{
    NodeWeight res;
    for(int i = 0; i < dimension; i++){
      res.weights[i] = weights[i] + ew.weights[i];
    }
    return res;
  }

  NodeWeight operator +=(const NodeWeight &ew){
    for(int i = 0; i < dimension; i++){
      weights[i] = weights[i] + ew.weights[i];
    }
    return *this;
  }

  bool operator!= (const NodeWeight &nw){
    for(int i = 0; i < dimension; i++){
      if(weights[i] != nw.weights[i]){
        return true;
      }
    }
    return false;
  }

  void operator -=(const NodeWeight &ew){
    for(int i = 0; i < dimension; i++){
      weights[i] = weights[i] - ew.weights[i];
    }
  }

  NodeWeight operator -(){
    NodeWeight res;
    for(int i = 0; i < dimension; i++){
      res.weights[i] = -weights[i];
    }
  }

  NodeWeight operator-(const NodeWeight &ew) const{
    NodeWeight res;
    for(int i = 0; i < dimension; i++){
        res.weights[i] = weights[i] - ew.weights[i];
    }
    return res;
  }

  NodeWeight operator*(const NodeWeight &ew) const{
    NodeWeight res;
    for(int i = 0; i < dimension; i++){
      res.weights[i] = weights[i] * ew.weights[i];
    }
    return res;
  }

  NodeWeight operator*(double d) const{
    NodeWeight res;
    for(int i = 0; i < dimension; i++){
      res.weights[i] = std::ceil(weights[i] * d);
    }
    return res;
  }

  NodeWeight operator*(std::array<double, dimension> d) const{
    NodeWeight res;
    for(int i = 0; i < dimension; i++){
      res.weights[i] = std::ceil(weights[i] * d[i]);
    }
    return res;
  }

  NodeWeight operator / (const NodeWeight &ew){
    NodeWeight res;
    for(int i = 0; i < dimension; i++){
      res.weights[i] = weights[i] / ew.weights[i];
    }
    return res;
  }

  bool operator<=(const NodeWeight &ew) const {
    for(int i = 0; i < dimension; i++){
      if(weights[i] > ew.weights[i]){
        return false;
      }
    }
    return true;
  }

  bool operator>=(const NodeWeight &ew) const {
    for(int i = 0; i < dimension; i++){
      if(weights[i] < ew.weights[i]){
        return false;
      }
    }
    return true;
  }

  bool operator > (int cmp) const { 
    for(int i = 0; i < dimension; i++){
      if(weights[i] > cmp){
        return true;
      }
    }
    return false;
  }

  bool operator == (int cmp) const {
    for(int i = 0; i < dimension; i++){
      if(weights[i] != cmp){
        return false;
      }
    }
    return true;
  }
  std::array<double,dimension> divide_with_double(double d){
    std::array<double,dimension> res;
    for(int i = 0; i < dimension; i++){
      res[i] = static_cast<double>(weights[i] / d);
    }
    return res;
  }

  NodeWeight min(const NodeWeight &nw) const{
    NodeWeight res;
    for(int i = 0; i < dimension; i++){
      res.weights[i] = std::min(weights[i], nw.weights[i]);
    }
    return res;
  }

  int sum() const{
    int sum = 0;
    for(int i = 0; i < dimension; i++){
      sum += weights[i];
    }
    return sum;
  }

  void store(int val){
    for(int i = 0; i < dimension; i++){
        weights[i] = 0;
    }
      
  }

  bool existsSmaller(int cmp){
    for(int i = 0; i < dimension; i++){
      if(weights[i] <= cmp){
        return true;
      }
    }
    return false;
  }


  

  NodeWeight add_fetch(const NodeWeight &ew, std::memory_order order){
    NodeWeight res;
    for(int i = 0; i < dimension; i++){
      res.weights[i] = __atomic_add_fetch(&weights[i], ew.weights[i], order);
      
    }
    return res;
  }

  NodeWeight fetch_add(const NodeWeight &ew, std::memory_order order){
    NodeWeight res;
    for(int i = 0; i < dimension; i++){
      res.weights[i] = __atomic_fetch_add(&weights[i], ew.weights[i], order);
      
    }
    return res;
  }

  NodeWeight fetch_sub(const NodeWeight &ew, std::memory_order order){
    NodeWeight res;
    for(int i = 0; i < dimension; i++){
      res.weights[i] = __atomic_fetch_sub(&weights[i], ew.weights[i], order);
      
    }
    return res;
  }

  NodeWeight sub_fetch(const NodeWeight &ew, std::memory_order order){
    NodeWeight res;
    for(int i = 0; i < dimension; i++){
      res.weights[i] = __atomic_sub_fetch(&weights[i], ew.weights[i], order);
      
    }
    return res;
  }

  double* multiply(double d[]){
    double res[dimension];
    for(int i = 0; i < dimension; i++){
      res[i] = static_cast<double>(weights[i]) * d[i];
    }
    return res;
  }

  uint32_t totalWeight() const{
    uint32_t res = 0;
    for(int i = 0; i < dimension; i++){
      res += weights[i];
    }
    return res;
  }

  double* norm(){
    double norm[dimension];
    uint32_t weight = totalWeight();
    for(int i = 0; i < dimension; i++){
      norm[i] = weights[i] /weight;
    }
    return norm;
  }  

  NodeWeight multiply(double d){
    NodeWeight res;
    for(int i = 0; i < dimension; i++){
      res.weights[i] = std::floor(d * weights[i]);
    }
  }

  bool operator>(std::array<double, dimension> d){
    for(int i = 0; i < dimension; i++){
      if(weights[i] <= d[i]){
        return false;
      }
    }
    return true;
  }
  
  bool chooseMoreBalanced(const NodeWeight &w1, const NodeWeight &w2, const NodeWeight &limit) const{
    int max_w1 = 0;
    int overflow_w1 = 0;
    int max_w2 = 0;
    int overflow_w2 = 0;
    for(int i = 0; i < dimension; i++){
      int overflow = weights[i] * (w1.weights[i] + weights[i] - limit.weights[i]);
      if(overflow > 0){
        overflow_w1 += overflow;
        max_w1 = std::max(max_w1, overflow_w1);
      }
      overflow = weights[i] * (w2.weights[i] + weights[i] - limit.weights[i]);
      if(overflow > 0){
        overflow_w2 += overflow;
        max_w2 = std::max(max_w2, overflow_w1);
      }
    }
    if(max_w1 != max_w2){
      return max_w1 < max_w2;
    }
    return overflow_w1 < overflow_w2;

  }

  bool isLighterPartition(const NodeWeight &p1, const NodeWeight &p2, const NodeWeight &max1, const NodeWeight &max2) const{
    double w1 = 0;
    double w2 = 0;
    for(int i = 0; i < dimension; i++){
      w1 += (max1.weights[i] - p1.weights[i]) / max1.weights[i];
      w2 += (max2.weights[i] - p2.weights[i]) / max2.weights[i];
    }
    return w1 > w2;
  }

  bool operator==(const NodeWeight &nw) const{
    for(int i = 0; i < dimension; i++){
      if(weights[i] != nw.weights[i]){
        return false;
      }
    }
    return true;
  }

  std::array<double,dimension> operator/(const std::array<double,dimension> d)const{
    std::array<double,dimension> res;
    for(int i = 0; i < dimension; i++){
      res[i] = weights[i] / d[i];
    }
    return res;
  }

  std::array<double, dimension> operator/(const NodeWeight &nw) const{
    std::array<double, dimension> res;
    for(int i = 0; i < dimension; i++){
      res[i] = weights[i] / static_cast<double>(nw.weights[i]);
    }
  }

  NodeWeight div(size_t div) const{
    NodeWeight res;
    for(int i = 0; i < dimension; i++){
      res.weights[i] = std::floor(weights[i] / div);
    }
    return res;
  }

  std::string to_string()const{
    std::string s;
    for(int i = 0; i < dimension; i++){
      s += weights[i];
      s += ' ';
    }
    return s;
  }

  double maxValues(const NodeWeight &nw) const{
    int d = 0;
    for(int i = 0; i < dimension; i++){
      d += std::max(weights[i], nw.weights[i]);
    }
    return d;    
  }

  void operator=(std::array<double, dimension> d){
    for(int i = 0; i < dimension; i++){
      weights[i] = std::ceil(d[i]);
    }
  }

  int weight() const{
    int res = 0;
    for(int i = 0; i < dimension; i++){
      res += weights[i];
    }
    return res;
  }

};

uint32_t scalar(NodeWeight &w1, NodeWeight &w2);

std::array<double, dimension> operator*(const double d, const NodeWeight &w);

std::ostream& operator<<(std::ostream& os,const  NodeWeight &nw);

void operator<<(std::ostringstream os, const std::array<double, dimension> arr);

bool equals_in_one_dimension(const NodeWeight &w1, const NodeWeight &w2);

std::array<double, dimension> operator+(double d, std::array<double, dimension> nw);

std::array<double, dimension> operator*(std::array<double, dimension> d, const NodeWeight &nw);

bool operator<=(const std::array<double, dimension> d1, const double d2[dimension]);

std::array<double, dimension> operator-(std::array<double, dimension> d, const NodeWeight &nw);

std::array<double, dimension> divide_to_double(const NodeWeight &nw1, const NodeWeight &nw2);

std::string to_string(const std::array<double, dimension> d);

using HardwareTopology = mt_kahypar::parallel::HardwareTopology<>;
using TBBInitializer = mt_kahypar::parallel::TBBInitializer<HardwareTopology, false>;

#define UI64(X) static_cast<uint64_t>(X)

#ifndef PARSER_H
#define PARSER_H
struct parallel_tag_t { };
#endif
using RatingType = double;
#if KAHYPAR_USE_64_BIT_IDS
#define ID(X) static_cast<uint64_t>(X)
using HypernodeID = uint64_t;
using HyperedgeID = uint64_t;
#else
#define ID(X) static_cast<uint32_t>(X)
using HypernodeID = uint32_t;
using HyperedgeID = uint32_t;
#endif
using HypernodeWeight = NodeWeight;
using HyperedgeWeight = int32_t;
using PartitionID = int32_t;
using Gain = int32_t;

// Graph Types
using NodeID = uint32_t;
using ArcWeight = double;

struct Arc {
  NodeID head;
  ArcWeight weight;

  Arc() :
    head(0),
    weight(0) { }

  Arc(NodeID head, ArcWeight weight) :
    head(head),
    weight(weight) { }
};

// Constant Declarations

static constexpr PartitionID kInvalidPartition = -1;
static constexpr HypernodeID kInvalidHypernode = std::numeric_limits<HypernodeID>::max();
static constexpr HypernodeID kInvalidHyperedge = std::numeric_limits<HyperedgeID>::max();
static constexpr Gain kInvalidGain = std::numeric_limits<HyperedgeID>::min();
static constexpr size_t kEdgeHashSeed = 42;

static constexpr HypernodeID invalidNode = std::numeric_limits<HypernodeID>::max();
static constexpr Gain invalidGain = std::numeric_limits<Gain>::min();

namespace ds {
  using Clustering = vec<PartitionID>;
  using Base=std::__atomic_base<uint32_t>;
}

struct Move {
  PartitionID from = kInvalidPartition;
  PartitionID to = kInvalidPartition;
  HypernodeID node = invalidNode;
  Gain gain = invalidGain;

  bool isValid() const {
    return from != kInvalidPartition;
  }

  void invalidate() {
    from = kInvalidPartition;
  }
};

struct Move_with_transformed_gain {
  PartitionID from = kInvalidPartition;
  PartitionID to = kInvalidPartition;
  HypernodeID node = invalidNode;
  double gain = static_cast<double>(invalidGain);

  bool isValid() const {
    return from != kInvalidPartition;
  }

  void invalidate() {
    from = kInvalidPartition;
  }
};

  struct Move_Internal{
    double gain_and_balance;
    Gain gain;
    double balance;

    bool operator<(const Move_Internal move) const{
      return gain_and_balance < move.gain_and_balance || gain_and_balance == move.gain_and_balance && gain > move.gain;
    }
    bool operator<=(const Move_Internal move) const{
      return gain_and_balance < move.gain_and_balance || gain_and_balance == move.gain_and_balance && gain >= move.gain;
    }
    bool operator>(const Move_Internal move) const{
      return gain_and_balance > move.gain_and_balance;
    }

    bool operator==(const Move_Internal move) const{
      return gain_and_balance == move.gain_and_balance && gain == move.gain && balance == move.balance;
    }

    void recomputeBalance(){
    if(balance != 0.0){
      gain_and_balance = gain > 0 ? gain * balance : gain / balance;
    }
    else{
      gain_and_balance = gain > 0 ? 0.0 : std::numeric_limits<double>::max();
    }
    }

    void setGain(Gain g){
      gain = g;
    }

    void addToGain(Gain g){
      __atomic_add_fetch(&gain, g, std::memory_order_relaxed);
    }
  };

  struct Rebalance_Move : Move_Internal{
    bool is_positive_move(){
      return balance < 0.0 || balance <= 0.0 && gain > 0;
    }
    
  };

  struct Refinement_Move : Move_Internal{
    
    bool is_positive_move(){
      return gain > 0 || gain == 0 && balance <= 0.0;
    }
  };

  /*struct FallBackMove{
    double gain_and_balance;
    Gain gain;
    double balance;

    bool operator<(const FallBackMove move) const{
      return (balance <= 0) && (move.balance <= 0) && gain_and_balance < move.gain_and_balance || !((balance <= 0) && (move.balance <= 0)) && balance < move.balance;
    }
    bool operator>(const FallBackMove move) const{
      return (balance <= 0) && (move.balance <= 0) && gain_and_balance > move.gain_and_balance || !((balance <= 0) && (move.balance <= 0)) && balance > move.balance;
    }

    bool operator==(const FallBackMove move) const{
      return gain_and_balance == move.gain_and_balance && gain == move.gain && balance == move.balance;
    }

    void recomputeBalance(){
      Gain tmp_gain = gain;
      double tmp_balance = balance;
      if(gain <= 0){
        tmp_gain -= 1;
      }
      if(balance <= 0){
        tmp_balance -= 0.0001;
      }
      gain_and_balance = tmp_gain > 0 ? -tmp_gain / tmp_balance : -tmp_gain * tmp_balance;
    }

    bool is_positive_move(){
      return true;
    }

    void setGain(Gain g){
      gain = g;
    }

    void addToGain(Gain g){
      __atomic_add_fetch(&gain, g, std::memory_order_relaxed);
    }
  };*/

  /*template<typename id, typename data>
  struct AddressablePQ{
    std::vector<id> v;
    std::vector<id> references;
    std::vector<bool> in_use;
    std::vector<data> gains_and_balances;
    uint64_t swap_up_counter = 0;
    uint64_t swap_down_counter = 0;
    void setSize(id size){
      references.resize(size);
      gains_and_balances.resize(size);
      in_use.resize(size);
      for(id i = 0; i < size; i++){
        in_use[i] = false;
      }
    }
    id getParent(id index){
      return (index + 1)/ 2 - 1;
    }
    id getChild(id index){
      return 2 * index;
    }
    double priority(id index){
      ASSERT(index < v.size());
      ASSERT(v[index] < gains_and_balances.size());
      return -gains_and_balances[v[index]].gain_and_balance;
    }
    void insert(std::pair<id, data> x){
      ASSERT(x.first < references.size());
      ASSERT(x.first != -1);
      if(!in_use[x.first]){
        id index = v.size();
        references[x.first] = index;
        in_use[x.first] = true;
        gains_and_balances[x.first] = x.second;
        v.push_back(x.first);
        siftUp(index);
      }
      else{
        ASSERT(references.size() > x.first);
        gains_and_balances[references[x.first]] = x.second;
        check(x.first);
      }
    }
    void siftDown(id index){
      id child = getChild(index);
      while (child < v.size()){
        id new_idx = child + static_cast<id>((child + 1 < v.size()) && (in_use[v[child]] < in_use[v[child + 1]]) 
          && (gains_and_balances[v[child]] > gains_and_balances[v[child + 1]]));
        if(!(in_use[v[index]] < in_use[v[new_idx]] || gains_and_balances[v[index]] > gains_and_balances[v[new_idx]] 
          && in_use[v[new_idx]])){
            break;
          }
        swap_down_counter++;
        references[v[index]] = new_idx;
        references[v[new_idx]] = index;
        std::swap(v[index], v[new_idx]);
        index = new_idx;
        child = getChild(index);
      }
    }
    void siftUp(id index){
      id parent = getParent(index);
      while (index > 0 && in_use[v[index]] && gains_and_balances[v[parent]] > gains_and_balances[v[index]]){
        references[v[index]] = parent;
        references[v[parent]] = index;
        std::swap(v[index], v[parent]);
        index = parent;
        parent = getParent(index);
        swap_up_counter++;
      }
    }
    std::pair<id, data> getMax(){
      ASSERT(v.size() > 0);
      return {v[0], gains_and_balances[v[0]]};
    }
    std::pair<id, data> deleteMax(){
      ASSERT(v.size() > 0);
      std::pair<id, data> max = getMax();
      v[0] = v[v.size() - 1];
      references[v[0]] = 0;
      v.pop_back();      
      in_use[max.first] = false;
      siftDown(0);
      return max;
    }
    data& get(id x){
      ASSERT(x < gains_and_balances.size());
      return gains_and_balances[x];
    }
    void updateKey(id index){
      ASSERT(index < gains_and_balances.size());
      gains_and_balances[index].recomputeBalance();
      insert({index, gains_and_balances[index]});
    }
    bool isEmpty(){
      return v.size() == 0;
    }
    void disable(id index){
      if(in_use[index]){
        in_use[index] = false;
        v[references[index]] = v[v.size() - 1];
        references[v[references[index]]] = references[index];
        v.pop_back();
        siftDown(references[index]);
      }      
    }
    void enable(int index){
      if(!in_use[index]){
        in_use[index] = true;
        references[index] = v.size();
        v.push_back(index);      
        siftUp(references[index]);
      }
    }
    bool isEnabled(id index){
      ASSERT(index < in_use.size());
      return in_use[index];
    }
    void check(id index){
      if(references[index] > 0 && gains_and_balances[index] < gains_and_balances[getParent(references[index])]){
        siftUp(references[index]);
      }
      siftDown(references[index]);
      
    }
    void sort(){
      for(id x = 1; x < v.size(); x++){
        siftUp(x);
      }
    }
  };*/

  /*struct Queues_per_node{
    using Balance=double;
    AddressablePQ<HypernodeID, Balance> top_moves;
    std::vector<std::pair<Balance, AddressablePQ<PartitionID, std::pair<Gain, Balance>>>> queues_per_node;
    void setIsolatedBalanceGain(HypernodeID hn, Balance balance_gain){
      queues_per_node[hn].first = balance_gain;
    }
    void update(HypernodeID hn){

    }

  };*/
  /*template<typename Move>
  struct AbstractMoveQueue{
    using Gain_and_Balance=double;
    using  Balance=double;
    AddressablePQ<HypernodeID, Gain_and_Balance> top_moves;
    std::vector<AddressablePQ<PartitionID, Move>> queues_per_node;
    std::vector<bool> locked;
    void initialize(HypernodeID num_nodes, PartitionID num_parts){
      top_moves.setSize(num_nodes);
      queues_per_node.resize(num_nodes);
      locked.resize(num_nodes);
      for(HypernodeID i = 0; i < num_nodes; i++){
        queues_per_node[i].setSize(num_parts);
      }
    }
    std::pair<HypernodeID, std::pair<PartitionID, Move>> deleteMax(){
      ASSERT(top_moves.v.size() > 0);
      HypernodeID node = top_moves.deleteMax().first;
      std::pair<PartitionID, Move> move = queues_per_node[node].deleteMax();
      if(!queues_per_node[node].isEmpty()){
        top_moves.insert({node, queues_per_node[node].getMax().second.gain_and_balance});
      }      
      return{node, move};
    }

    std::pair<HypernodeID, std::pair<PartitionID, Move>> getMax(){
      ASSERT(top_moves.v.size() > 0);
      HypernodeID node = top_moves.getMax().first;
      std::pair<PartitionID, Move> move = queues_per_node[node].getMax();   
      return{node, move};
    }

    void insert(std::pair<HypernodeID, std::pair<PartitionID, Move>> move){
      if(move.second.second.is_positive_move()){
        ASSERT(queues_per_node.size() > move.first);
        queues_per_node[move.first].insert(move.second);
        update(move.first);
      }
      else{
        ASSERT(queues_per_node.size() > move.first);
        ASSERT(queues_per_node[move.first].gains_and_balances.size() > move.second.first);
        queues_per_node[move.first].gains_and_balances[move.second.first] = move.second.second;
        if(queues_per_node[move.first].isEnabled(move.second.first)){
          disable({move.first, move.second.first});
        }
      }
    }
    void insert_without_updating(std::pair<HypernodeID, std::pair<PartitionID, Move>> move, bool disabled = false){
      
      if(move.second.second.is_positive_move() && !disabled){
        ASSERT(queues_per_node.size() > move.first);
        queues_per_node[move.first].insert(move.second);
      }
      else{
        ASSERT(queues_per_node.size() > move.first);
        ASSERT(queues_per_node[move.first].gains_and_balances.size() > move.second.first);
        queues_per_node[move.first].gains_and_balances[move.second.first] = move.second.second;
        if(queues_per_node[move.first].isEnabled(move.second.first)){
          queues_per_node[move.first].disable(move.second.first);
        }
      }
    }
    bool isEmpty(){
      return top_moves.isEmpty();
    }

    Move get(std::pair<HypernodeID, PartitionID> move){
      return queues_per_node[move.first].get(move.second);
    }

    void lock(HypernodeID hn){
      locked[hn] = true;
    }

    
    void addToGain(std::pair<HypernodeID, std::pair<PartitionID, Gain>> x){
      queues_per_node[x.first].gains_and_balances[x.second.first].addToGain(x.second.second);
    }
    
    void update(HypernodeID hn){
      if(queues_per_node[hn].v.size() == 0 || locked[hn] == true){
        top_moves.disable(hn);
      }
      else if(!top_moves.isEnabled(hn) || top_moves.get(hn) != queues_per_node[hn].getMax().second.gain_and_balance){
        top_moves.insert({hn, queues_per_node[hn].getMax().second.gain_and_balance});
      }
    }
    void disable(std::pair<HypernodeID, PartitionID> x){
      queues_per_node[x.first].disable(x.second);
      update(x.first);
    }
    void checkSizes(PartitionID size){
      for(HypernodeID hn = 0; hn < queues_per_node.size(); hn++){
        if(queues_per_node[hn].v.size() <= size);
      }
    }
    void check(){
      for(HypernodeID hn = 0; hn < queues_per_node.size(); hn++){
        update(hn);
      }
    }
    void resetGainAndDisable(std::pair<HypernodeID, PartitionID> x){
        queues_per_node[x.first].gains_and_balances[x.second].gain = 0;
        queues_per_node[x.first].gains_and_balances[x.second].balance = 0.0;
        disable(x);
    }
  };*/

  /*struct MoveQueue : AbstractMoveQueue<Move_internal> {
    void changeBalance(std::pair<HypernodeID, std::pair<PartitionID, Balance>> x){
      queues_per_node[x.first].gains_and_balances[x.second.first].balance = x.second.second;
      if(!queues_per_node[x.first].gains_and_balances[x.second.first].is_positive_move()){
        if(queues_per_node[x.first].isEnabled(x.second.first)){
          queues_per_node[x.first].disable(x.second.first);
        }
      }
      else{
        queues_per_node[x.first].gains_and_balances[x.second.first].recomputeBalance();
        if(!queues_per_node[x.first].isEnabled(x.second.first)){
          queues_per_node[x.first].enable(x.second.first);
        }
        else{
          queues_per_node[x.first].check(x.second.first);
        }
      }
      
    }

    void changeGain(std::pair<HypernodeID, std::pair<PartitionID, Gain>> x){
      queues_per_node[x.first].gains_and_balances[x.second.first].setGain(x.second.second);
      if(!queues_per_node[x.first].gains_and_balances[x.second.first].is_positive_move()){
        if(queues_per_node[x.first].isEnabled(x.second.first)){
          queues_per_node[x.first].disable(x.second.first);
        }
      }
      else{
        queues_per_node[x.first].gains_and_balances[x.second.first].recomputeBalance();
        if(!queues_per_node[x.first].isEnabled(x.second.first)){
          queues_per_node[x.first].enable(x.second.first);
        }
        else{
          queues_per_node[x.first].check(x.second.first);
        }
      }
    }
  };*/ 

  /*struct Fallback_MoveQueue : AbstractMoveQueue<FallBackMove>{
    void changeBalance(std::pair<HypernodeID, std::pair<PartitionID, Balance>> x){
      queues_per_node[x.first].gains_and_balances[x.second.first].balance = x.second.second;
      queues_per_node[x.first].gains_and_balances[x.second.first].recomputeBalance();
        if(!queues_per_node[x.first].isEnabled(x.second.first)){
          queues_per_node[x.first].enable(x.second.first);
        }
        else{
          queues_per_node[x.first].check(x.second.first);
        }
    }
  };*/


struct Memento {
  HypernodeID u; // representative
  HypernodeID v; // contraction partner
};

template<typename Hypergraph>
struct ExtractedHypergraph {
  Hypergraph hg;
  vec<HypernodeID> hn_mapping;
  vec<uint8_t> already_cut;
};

using Batch = parallel::scalable_vector<Memento>;
using BatchVector = parallel::scalable_vector<Batch>;
using VersionedBatchVector = parallel::scalable_vector<BatchVector>;

using MoveID = uint32_t;
using SearchID = uint32_t;

// Forward Declaration
class TargetGraph;
namespace ds {
class Bitset;
class StaticGraph;
class PinCountSnapshot;
class StaticHypergraph;
class DynamicGraph;
class DynamicHypergraph;
class ConnectivityInfo;
class SparseConnectivityInfo;
}

struct SynchronizedEdgeUpdate {
  HyperedgeID he = kInvalidHyperedge;
  PartitionID from = kInvalidPartition;
  PartitionID to = kInvalidPartition;
  HyperedgeID edge_weight = 0;
  HypernodeID edge_size = 0;
  HypernodeID pin_count_in_from_part_after = kInvalidHypernode;
  HypernodeID pin_count_in_to_part_after = kInvalidHypernode;
  PartitionID block_of_other_node = kInvalidPartition;
  mutable ds::Bitset* connectivity_set_after = nullptr;
  mutable ds::PinCountSnapshot* pin_counts_after = nullptr;
  const TargetGraph* target_graph = nullptr;
  ds::Array<SpinLock>* edge_locks = nullptr;
};

struct NoOpDeltaFunc {
  void operator() (const SynchronizedEdgeUpdate&) { }
};

template<typename Hypergraph, typename ConInfo>
struct PartitionedHypergraphType {
  static constexpr mt_kahypar_partition_type_t TYPE = NULLPTR_PARTITION;
};

template<>
struct PartitionedHypergraphType<ds::StaticHypergraph, ds::ConnectivityInfo> {
  static constexpr mt_kahypar_partition_type_t TYPE = MULTILEVEL_HYPERGRAPH_PARTITIONING;
};

template<>
struct PartitionedHypergraphType<ds::StaticHypergraph, ds::SparseConnectivityInfo> {
  static constexpr mt_kahypar_partition_type_t TYPE = LARGE_K_PARTITIONING;
};

template<>
struct PartitionedHypergraphType<ds::DynamicHypergraph, ds::ConnectivityInfo> {
  static constexpr mt_kahypar_partition_type_t TYPE = N_LEVEL_HYPERGRAPH_PARTITIONING;
};

template<typename Graph>
struct PartitionedGraphType {
  static constexpr mt_kahypar_partition_type_t TYPE = NULLPTR_PARTITION;
};

template<>
struct PartitionedGraphType<ds::StaticGraph> {
  static constexpr mt_kahypar_partition_type_t TYPE = MULTILEVEL_GRAPH_PARTITIONING;
};

template<>
struct PartitionedGraphType<ds::DynamicGraph> {
  static constexpr mt_kahypar_partition_type_t TYPE = N_LEVEL_GRAPH_PARTITIONING;
};




  template<typename id, typename data>
  struct AddressablePQ{
    std::vector<id> pq;
    std::vector<id> references;
    std::vector<data> items;
    std::vector<bool> in_use;
    void setSize(id size){
      references.resize(size);
      items.resize(size);
      in_use.resize(size, false);
    }

    bool isInRange(id index){
      return index >= 0 && index < pq.size();
    }
    id getParent(id index){
      return (index - 1)/ 2;
    }
    id getChild(id index){
      return 2 * index + 1;
    }
    void insert(std::pair<id, data> x){
      references[x.first] = pq.size();
      pq.push_back(x.first);
      items[x.first] = x.second;
      in_use[x.first] = true;
      siftUp(pq.size() - 1);
    }

    void swapInPQ(size_t index1, size_t index2){
      ASSERT(references[pq[index1]] == index1);
      ASSERT(references[pq[index2]] == index2);
      references[pq[index1]] = index2;
      references[pq[index2]] = index1;
      size_t tmp = pq[index1];
      pq[index1] = pq[index2];
      pq[index2] = tmp;
    }
    
    void siftUp(size_t index){
      if(index == 0){
        return;
      }
      size_t parent = getParent(index);
      if(items[pq[index]] < items[pq[parent]]){
        swapInPQ(index, parent);
        siftUp(parent);
      }
    }

    void siftDown(size_t index){
      if(getChild(index) >= pq.size()) return;
      size_t child1 = getChild(index);
      size_t child2 = std::min(child1 + 1, pq.size() - 1);
      size_t chosen_child = child1 + static_cast<size_t>((items[pq[child1]] > items[pq[child2]]));
      if(items[pq[index]] > items[pq[chosen_child]]){
        swapInPQ(index, chosen_child);
        siftDown(chosen_child); 
      }
    }

    void changeKey(std::pair<id,data> pair){
      if(!in_use[pair.first]){
        insert(pair);
      }
      else{
        if(pair.second < items[pair.first]){
          items[pair.first] = pair.second;
          siftUp(references[pair.first]);
        }
        else if(pair.second > items[pair.first]){
          items[pair.first] = pair.second;
          siftDown(references[pair.first]);
        }
      }
      
    }

    std::pair<id, data> getMax(){
      return {pq[0], items[pq[0]]};
    }

    std::pair<id, data> deleteMax(){
      auto tmp = std::pair<id, data>(pq[0], items[pq[0]]);
      invalidate(pq[0]);
      return tmp;
    }

    void invalidate(id idx){
      if(in_use[idx]){
        in_use[idx] = false;
        pq[references[idx]] = pq[pq.size() - 1];
        references[pq[references[idx]]] = references[idx];
        pq.pop_back();
        if(references[idx] < pq.size()){
          siftDown(references[idx]);
        }
      }
    }

    bool isEmpty(){
      return pq.size() == 0;
    }

    data get(id idx){
      return items[idx];
    }

  };






































}
// namespace mt_kahypar
/*namespace std{
  ostringstream& operator<<(ostringstream& s, const array<double,mt_kahypar::dimension> arr){
    for(int i = 0; i < mt_kahypar::dimension; i++){
      s << arr[i] << ' ';
    }
    return s;
  }
}*/
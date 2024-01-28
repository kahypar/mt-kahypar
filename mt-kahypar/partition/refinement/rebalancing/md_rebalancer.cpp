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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/cast.h"

namespace mt_kahypar{

  

  struct Move_internal{
    
  double gain_and_balance;
  double gain;
  double balance;

    bool operator<(const Move_internal move) const{
      return gain_and_balance < move.gain_and_balance;
    }
    bool operator>(const Move_internal move) const{
      return gain_and_balance > move.gain_and_balance;
    }

    bool operator==(const Move_internal move) const{
      return gain_and_balance == move.gain_and_balance && gain == move.gain && balance == move.balance;
    }

    void recomputeBalance(){
    Gain tmp_gain = gain;
    if(gain <= 0){
      tmp_gain -= 0.01;
    }
    gain_and_balance = tmp_gain > 0 ? -tmp_gain / (balance - 0.001) : -tmp_gain * (balance - 0.001);
    }
  };

  template<typename id, typename data>
  struct AddressablePQ{
    std::vector<id> v;
    std::vector<id> references;
    std::vector<bool> in_use;
    std::vector<data> gains_and_balances;
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
        if(references[x.first] > 0 && x.second < gains_and_balances[getParent(references[x.first])]){
          siftUp(references[x.first]);
        }
        else if(getChild(references[x.first]) < v.size() && x.second > gains_and_balances[getChild(references[x.first])]){
          siftDown(references[x.first]);
        }
      }
    }
    void siftDown(id index){
      id child = getChild(index);
      while (child < v.size() && (gains_and_balances[v[child]] < gains_and_balances[v[index]] || !in_use[v[child]])){
        id new_idx = child + static_cast<id>((child + 1 < v.size()) ? 
          (gains_and_balances[v[child]] > gains_and_balances[v[child + 1]]) : 0);
        references[v[index]] = new_idx;
        references[v[new_idx]] = index;
        std::swap(v[index], v[new_idx]);
        index = new_idx;
        child = getChild(index);
      }
    }
    void siftUp(id index){
      id parent = getParent(index);
      while (index > 0 && gains_and_balances[v[parent]] > gains_and_balances[v[index]]){
        references[v[index]] = parent;
        references[v[parent]] = index;
        std::swap(v[index], v[parent]);
        index = parent;
        parent = getParent(index);
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
      in_use[index] = false;
      siftDown(index);
    }
    bool isEnabled(id index){
      ASSERT(index < in_use.size());
      return in_use[index];
    }
  };

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
  struct MoveQueue{
    using Gain_and_Balance=double;
    AddressablePQ<HypernodeID, Gain_and_Balance> top_moves;
    std::vector<AddressablePQ<PartitionID, Move_internal>> queues_per_node;
    void initialize(HypernodeID num_nodes, PartitionID num_parts){
      top_moves.setSize(num_nodes);
      queues_per_node.resize(num_nodes);
      for(HypernodeID i = 0; i < num_nodes; i++){
        queues_per_node[i].setSize(num_parts);
      }
    }
    std::pair<HypernodeID, PartitionID> deleteMax(){
      HypernodeID node = top_moves.deleteMax().first;
      PartitionID to = queues_per_node[node].deleteMax().first;
      if(!queues_per_node[node].isEmpty()){
        top_moves.insert({node, queues_per_node[node].getMax().second.gain_and_balance});
      }      
      return{node, to};
    }
    void insert(std::pair<HypernodeID, std::pair<PartitionID, Move_internal>> move){
      if(move.second.second.balance <= 0){
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
    bool isEmpty(){
      return top_moves.isEmpty();
    }
    void changeBalance(std::pair<HypernodeID, std::pair<PartitionID, double>> x){
      Move_internal move = queues_per_node[x.first].get(x.second.first);
      move.balance = x.second.second;
      move.recomputeBalance();
      insert({x.first, {x.second.first, move}});
    }
    void addToGain(std::pair<HypernodeID, std::pair<PartitionID, Gain>> x){
      Move_internal move = queues_per_node[x.first].get(x.second.first);
      move.gain += x.second.second;
      move.recomputeBalance();
      insert({x.first, {x.second.first, move}});
    }
    void update(HypernodeID hn){
      ASSERT(queues_per_node[hn].v.size() > 0);
      if(top_moves.get(hn) != queues_per_node[hn].getMax().second.gain_and_balance){
        top_moves.insert({hn, queues_per_node[hn].getMax().second.gain_and_balance});
      }
    }
    void disable(std::pair<HypernodeID, PartitionID> x){
      queues_per_node[x.first].disable(x.second);
      update(x.first);
    }
  };

  template <typename GraphAndGainTypes>
  bool MDRebalancer<GraphAndGainTypes>::refineInternal(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                       vec<vec<Move>>* moves_by_part,
                                                       vec<Move>* moves_linear,
                                                       Metrics& best_metrics) {
                                                      
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
    }

    auto balance_gain = [&](const PartitionedHypergraph& phg, HypernodeID node, PartitionID from, PartitionID to){
      double gain = 0.0;
      for(int i = 0; i < dimension; i++){
        gain += std::max(0, std::min(phg.nodeWeight(node).weights[i], phg.partWeight(to).weights[i] + 
        phg.nodeWeight(node).weights[i] - _context.partition.max_part_weights[to].weights[i])) * _context.partition.max_part_weights_inv[to][i]
        - std::max(0, std::min(phg.nodeWeight(node).weights[i], phg.partWeight(from).weights[i] - _context.partition.max_part_weights[from].weights[i])) * _context.partition.max_part_weights_inv[from][i];
      }
      return gain;
    };
    MoveQueue queue;
    queue.initialize(phg.initialNumNodes(), phg.k());
    for(HypernodeID hn : phg.nodes()){
      std::vector<Move_md> moves = _gain.allGains(phg, hn);
      for(int i = 0; i < moves.size(); i++){
        ASSERT(phg.partID(hn) != moves[i].to);
        queue.insert({hn, {moves[i].to, {moves[i].gain_and_balance, moves[i].gain, moves[i].balance}}});
      }   
    }
    PartitionID imbalanced = 0;
    for(PartitionID p = 0; p < phg.k(); p++){
      if(phg.partWeight(p) > _context.partition.max_part_weights[p]){
        imbalanced++;
      }
    }
    while(imbalanced != 0 && !queue.isEmpty()){
      std::cout << "marker0\n";
      /*for(int i = 0; i < 1; i++){
        std::cout << queue.queues_per_node[queue.top_moves.v[0]].v[i] << "\n";
        std::cout << queue.queues_per_node[queue.top_moves.v[0]].gains_and_balances[i].gain_and_balance << "\n";
         std::cout << queue.queues_per_node[queue.top_moves.v[0]].gains_and_balances[i].gain << "\n";
          std::cout << queue.queues_per_node[queue.top_moves.v[0]].gains_and_balances[i].balance << "\n";
      }*/
      /*std::cout << phg.partID(queue.top_moves.v[0]) << "\n";
      for(int i = 0; i < queue.queues_per_node[queue.top_moves.v[0]].v.size(); i++){
        std::cout << queue.queues_per_node[queue.top_moves.v[0]].v[i] << "\n";
      }
      std::cout << "end1\n";
      for(int i = 0; i < queue.queues_per_node[queue.top_moves.v[0]].gains_and_balances.size(); i++){
        std::cout << queue.queues_per_node[queue.top_moves.v[0]].gains_and_balances[i].gain << "\n" 
        << queue.queues_per_node[queue.top_moves.v[0]].gains_and_balances[i].balance << "\n"
        << queue.queues_per_node[queue.top_moves.v[0]].gains_and_balances[i].gain_and_balance << "\n";
        std::cout << "\n";
      }*/
      std::pair<HypernodeID, PartitionID> max_move = queue.deleteMax();
      HypernodeID node = max_move.first;
      Move move = {phg.partID(node), max_move.second, node, 0};
      imbalanced += (phg.partWeight(phg.partID(node)) - phg.nodeWeight(node) > _context.partition.max_part_weights[phg.partID(node)])
        - (phg.partWeight(phg.partID(node)) > _context.partition.max_part_weights[phg.partID(node)])
        + (phg.partWeight(move.to) > _context.partition.max_part_weights[move.to])
        - (phg.partWeight(move.to) + phg.nodeWeight(node) > _context.partition.max_part_weights[move.to]);
      phg.changeNodePart(node, phg.partID(node), move.to, objective_delta);
      queue.disable({node, move.to});
      std::cout << "marker1\n";
      for(HypernodeID hn : phg.nodes()){
        std::cout << hn << "\n";
        if(phg.partID(hn) != move.to){
          queue.changeBalance({hn, {move.to, balance_gain(phg, hn, phg.partID(hn), move.to)}});
        }
        if(phg.partID(hn) != move.from){
          queue.changeBalance({hn, {move.from, balance_gain(phg, hn, phg.partID(hn), move.from)}});
        }
        if(phg.partID(hn) == move.from || phg.partID(hn) == move.to){
          for(PartitionID p = 0; p < phg.k(); p++){
            if(p != phg.partID(hn)){
              queue.changeBalance({hn, {p, balance_gain(phg, hn, phg.partID(hn), p)}});
            }
          }
        }
      }
      std::cout << "marker2\n";
      for(Move m : _gain.getChangedMoves(phg, {move.from, move.to, move.node, 0})){
        if(m.to != phg.partID(m.node)){
          queue.addToGain({m.node, {m.to, m.gain}});
        }        
      }
    }
    /*
    std::array<std::priority_queue<std::pair<int64_t, PartitionID>>, mt_kahypar::dimension> max_imbalances;
    for(int j = 0; j < mt_kahypar::dimension; j++){
      for(int i = 0; i < phg.k(); i++){    
        int64_t imbalance = phg.partWeight(i).weights[j] - _context.partition.max_part_weights[i].weights[j];
        if(imbalance > 0){
          max_imbalances[j].push(std::pair<int64_t, PartitionID>(imbalance, i));
        }        
      }
    }
    _gain.reset();
    for(HypernodeID hn : phg.nodes()){
      Move move = _gain.basicMaxGainMove_global_gain(phg, hn, max_imbalances);
      if(move.to != move.from){
        phg.changeNodePart(hn, move.from, move.to, objective_delta);
      }
      for(int j = 0; j < mt_kahypar::dimension; j++){
        int64_t imbalance = phg.partWeight(move.to).weights[j] - _context.partition.max_part_weights[move.to].weights[j];
        if(imbalance > 0){
          max_imbalances[j].push(std::pair<int64_t, PartitionID>(imbalance, move.to));
        }
        
      }
      
    }*/

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

    bool improvement = delta < 0;
    return improvement;
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

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
    phg.doParallelForAllNodes([&](const HypernodeID hn){
      std::vector<Move_md> moves = _gain.allGains(phg, hn);
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
    int32_t counter = 0;
    while(imbalanced != 0 && !queue.isEmpty()){
      if(counter % 500 == 0){
        double ib = 0.0;
        for(int i = 0; i < phg.k(); i++){
          for(int j = 0; j < dimension; j++){
            ib += std::max(0.0, (phg.partWeight(i).weights[j] - _context.partition.max_part_weights[i].weights[j]) * _context.partition.max_part_weights_inv[i][j]);
          }
        }
        std::cout << ib << "\n";
      }
      counter++;  /*for(int i = 0; i < 1; i++){
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
      std::pair<HypernodeID, std::pair<PartitionID, Move_internal>> max_move = queue.deleteMax();
      ASSERT(phg.partID(max_move.first) != -1);
      ASSERT(max_move.second != -1);
      HypernodeID node = max_move.first;
      PartitionID to = max_move.second.first;
      Move_internal gains = max_move.second.second;
      Move move = {phg.partID(node), to, node, 0};
      imbalanced += (phg.partWeight(phg.partID(node)) - phg.nodeWeight(node) > _context.partition.max_part_weights[phg.partID(node)])
        - (phg.partWeight(phg.partID(node)) > _context.partition.max_part_weights[phg.partID(node)])
        - (phg.partWeight(to) > _context.partition.max_part_weights[to])
        + (phg.partWeight(to) + phg.nodeWeight(node) > _context.partition.max_part_weights[to]);
      queue.insert_without_updating({node, {phg.partID(node), {-gains.gain_and_balance, -gains.gain, -gains.balance}}});
      phg.changeNodePart(node, phg.partID(node), to, objective_delta);
      queue.resetGainAndDisable({node, to});
      _gain.getChangedMoves(phg, {move.from, move.to, move.node, 0}, &queue);
      /*for(Move m : _gain.getChangedMoves(phg, {move.from, move.to, move.node, 0})){
        if(m.to != phg.partID(m.node)){
          queue.addToGain({m.node, {m.to, m.gain}});
        }        
      }*/
      phg.doParallelForAllNodes([&](HypernodeID hn){
        if(phg.partID(hn) == move.from || phg.partID(hn) == move.to){
          for(PartitionID p = 0; p < phg.k(); p++){
            if(p != phg.partID(hn)){
              queue.changeBalance({hn, {p, balance_gain(phg, hn, phg.partID(hn), p)}});
            }
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
      });
      queue.check();
      queue.checkSizes(phg.k());
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

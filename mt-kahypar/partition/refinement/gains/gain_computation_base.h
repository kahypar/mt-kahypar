/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
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

#include <vector>

#include "kahypar-resources/meta/mandatory.h"

#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {

struct Move_md{
  double gain_and_balance;
  PartitionID to;
  double gain;
  double balance;

  bool operator<(const Move_md move) const{
    return gain_and_balance < move.gain_and_balance;
  }

  bool operator>(const Move_md move) const{
    return gain_and_balance > move.gain_and_balance;
  }

  bool operator==(const Move_md move) const{
      return gain_and_balance == move.gain_and_balance && gain == move.gain && balance == move.balance && to == move.to;
    }

  void recomputeBalance(){
    gain_and_balance = gain > 0 ? -gain / balance : -gain * balance;
  }
};

template <class Derived = Mandatory,
          class AttributedGains = Mandatory>
class GainComputationBase {
  using DeltaGain = tbb::enumerable_thread_specific<Gain>;

 public:
  using RatingMap = ds::SparseMap<PartitionID, Gain>;
  using TmpScores = tbb::enumerable_thread_specific<RatingMap>;

  GainComputationBase(const Context& context,
                      const bool disable_randomization) :
    _context(context),
    _disable_randomization(disable_randomization),
    _deltas(0),
    _tmp_scores([&] {
      return constructLocalTmpScores();
    }) { }

  template<typename PartitionedHypergraph>
  std::vector<Move_md> allGains(const PartitionedHypergraph& phg,
                          const HypernodeID hn) {
    Derived* derived = static_cast<Derived*>(this);
    RatingMap& tmp_scores = _tmp_scores.local();
    Gain isolated_block_gain = 0;
    tmp_scores.clear();
    std::vector<Move_md> gains;
    derived->precomputeGains(phg, hn, tmp_scores, isolated_block_gain, true);
    auto balance_gain = [&](const PartitionedHypergraph& phg, HypernodeID node, PartitionID from, PartitionID to){
      double gain = 0.0;
      for(int i = 0; i < dimension; i++){
        gain += std::max(0, std::min(phg.nodeWeight(node).weights[i], phg.partWeight(to).weights[i] + 
        phg.nodeWeight(node).weights[i] - _context.partition.max_part_weights[to].weights[i])) * _context.partition.max_part_weights_inv[to][i]
        - std::max(0, std::min(phg.nodeWeight(node).weights[i], phg.partWeight(from).weights[i] - _context.partition.max_part_weights[from].weights[i])) * _context.partition.max_part_weights_inv[from][i];
      }
      return gain;
    };
    /*auto balance_metric_gain = [&](const PartitionedHypergraph& phg, HypernodeID node, PartitionID from, PartitionID to){
      double gain = 0.0;
      for(int i = 0; i < dimension; i++){
        if(max_imbalance[])

      }
    }*/
    PartitionID from = phg.partID(hn);
    for(PartitionID p = 0; p < phg.k(); p++){
      if(p != from){
      double balance_new = balance_gain(phg, hn, from, p);
      double gain_new = tmp_scores[p];
      double gainandbalance = gain_new > 0 ? -gain_new / (balance_new - 0.001) : -gain_new * (balance_new - 0.001);
      gains.push_back({gainandbalance, p, gain_new, balance_new});
      }
    }
    return gains;

  }


  /*template<typename PartitionedHypergraph>
  Move_md basicMaxGainMove(const PartitionedHypergraph& phg,
                          const HypernodeID hn) {
    Derived* derived = static_cast<Derived*>(this);
    RatingMap& tmp_scores = _tmp_scores.local();
    Gain isolated_block_gain = 0;
    tmp_scores.clear();
    derived->precomputeGains(phg, hn, tmp_scores, isolated_block_gain, true);
    auto balance_gain = [&](const PartitionedHypergraph& phg, HypernodeID node, PartitionID from, PartitionID to){
      double gain = 0.0;
      for(int i = 0; i < dimension; i++){
        gain += std::max(0, std::min(phg.nodeWeight(node).weights[i], phg.partWeight(to).weights[i] + 
        phg.nodeWeight(node).weights[i] - _context.partition.max_part_weights[to].weights[i])) * _context.partition.max_part_weights_inv[to][i]
        - std::max(0, std::min(phg.nodeWeight(node).weights[i], phg.partWeight(from).weights[i] - _context.partition.max_part_weights[from].weights[i])) * _context.partition.max_part_weights_inv[from][i];
      }
      return gain;
    };
    auto balance_metric_gain = [&](const PartitionedHypergraph& phg, HypernodeID node, PartitionID from, PartitionID to){
      double gain = 0.0;
      for(int i = 0; i < dimension;ö, std::numeric_limits<double>().max()};
    for(PartitionID p = 0; p < phg.k() && p != from; p++){
      double balance_new = balance_gain(phg, hn, from, p);
      double gain_new = tmp_scores[p];
      double gainandbalance = gain_new > 0 ? -gain_new / balance_new : -gain_new * balance_new;
      if(balance_new < 0 && gainandbalance < best_move.first){
        best_move = {gainandbalance, hn, p, gain_new, balance_new};
      }
    }
    return best_move;

  }*/

  template<typename PartitionedHypergraph>
  void getChangedMoves(const PartitionedHypergraph& phg, const Move move, MoveQueue* mq){
    Derived* derived = static_cast<Derived*>(this);
    derived->getChangedMoves(phg, move, mq);
  }

  template<typename PartitionedHypergraph>
  Move basicMaxGainMove_global_gain(const PartitionedHypergraph& phg,
                          const HypernodeID hn, std::array<std::priority_queue<std::pair<int64_t, PartitionID>>, mt_kahypar::dimension> max_imbalances) {
    Derived* derived = static_cast<Derived*>(this);
    RatingMap& tmp_scores = _tmp_scores.local();
    tmp_scores.clear();
    PartitionID from = phg.partID(hn);
    double isolated_block_balance_gain = 0.0;
    auto pop_expired = [&](std::priority_queue<std::pair<int64_t, PartitionID>> pq, int dimension){
      while(pq.top().first != phg.partWeight(pq.top().second).weights[dimension] - _context.partition.max_part_weights[pq.top().second].weights[dimension]){
        pq.pop();
      }
    };
    for(int i = 0; i < mt_kahypar::dimension; i++){
      pop_expired(max_imbalances[i], i);
      if(max_imbalances[i].top().second == from){
        std::pair<int64_t, PartitionID> old_max = max_imbalances[i].top();
        max_imbalances[i].pop();
        int64_t balance_gain = std::min(static_cast<int64_t>(phg.nodeWeight(hn).weights[i]), old_max.first - (max_imbalances[i].size() > 0 ? max_imbalances[i].top().first : 0));
        isolated_block_balance_gain += balance_gain * _context.partition.max_part_weights_inv[from][i];
        int64_t new_imbalance = old_max.first - phg.nodeWeight(hn).weights[i];
        if(new_imbalance > 0){
          max_imbalances[i].push(std::pair<int64_t,PartitionID>(new_imbalance, from));
        }
        
      }
    }
    Gain isolated_block_gain = 0;
    derived->precomputeGains(phg, hn, tmp_scores, isolated_block_gain, true);
    auto balance_loss = [&](const PartitionedHypergraph& phg, HypernodeID node, PartitionID from, PartitionID to){
      double gain = 0.0;
      for(int i = 0; i < dimension; i++){
        gain += (phg.nodeWeight(node).weights[i] + phg.partWeight(to).weights[i] - _context.partition.max_part_weights[to].weights[i] - 
        max_imbalances[i].size() > 0 ? max_imbalances[i].top().first : 0) * _context.partition.max_part_weights_inv[to][i];
        
      }
      return gain;
    };
    /*auto balance_metric_gain = [&](const PartitionedHypergraph& phg, HypernodeID node, PartitionID from, PartitionID to){
      double gain = 0.0;
      for(int i = 0; i < dimension; i++){
        if(max_imbalance[])

      }
    }*/
    
    Move best_move { from, from, hn, 0};
    double gainandbalance_old = std::numeric_limits<double>().max();
    for(PartitionID p = 0; p < phg.k() && p != from; p++){
      double loss = balance_loss(phg, hn, phg.partID(hn), p);
      double gain_new = tmp_scores[p];
      double gainandbalance = gain_new > 0 ? -gain_new / (loss - isolated_block_gain) : -gain_new * (loss - isolated_block_gain);
      if(loss < isolated_block_gain  && gainandbalance < gainandbalance_old){
        gainandbalance_old = gainandbalance;
        best_move.to = p;
        best_move.gain = gain_new;
      }
    }
    return best_move;

  }






  template<typename PartitionedHypergraph>
  Move computeMaxGainMove(const PartitionedHypergraph& phg,
                          const HypernodeID hn,
                          const bool rebalance = false,
                          const bool consider_non_adjacent_blocks = false,
                          const bool allow_imbalance = false) {
    Derived* derived = static_cast<Derived*>(this);
    RatingMap& tmp_scores = _tmp_scores.local();
    Gain isolated_block_gain = 0;
    derived->precomputeGains(phg, hn, tmp_scores, isolated_block_gain, consider_non_adjacent_blocks);

    PartitionID from = phg.partID(hn);
    Move best_move { from, from, hn, rebalance ? std::numeric_limits<Gain>::max() : 0 };
    HypernodeWeight hn_weight = phg.nodeWeight(hn);
    int cpu_id = THREAD_ID;
    utils::Randomize& rand = utils::Randomize::instance();
    auto test_and_apply = [&](const PartitionID to,
                              const Gain score,
                              const bool no_tie_breaking = false) {
      bool new_best_gain = (score < best_move.gain) ||
                            (score == best_move.gain &&
                            !_disable_randomization &&
                            (no_tie_breaking || rand.flipCoin(cpu_id)));
      if (new_best_gain && (allow_imbalance || phg.partWeight(to) + hn_weight <=
          _context.partition.max_part_weights[to])) {
        best_move.to = to;
        best_move.gain = score;
        return true;
      } else {
        return false;
      }
    };

    for ( const auto& entry : tmp_scores ) {
      const PartitionID to = entry.key;
      if (from != to) {
        const Gain score = derived->gain(entry.value, isolated_block_gain);
        test_and_apply(to, score);
      }
    }

    if ( consider_non_adjacent_blocks && best_move.to == from ) {
      // This is important for our rebalancer as the last fallback strategy
      vec<PartitionID> non_adjacent_block;
      for ( PartitionID to = 0; to < _context.partition.k; ++to ) {
        if ( from != to && !tmp_scores.contains(to) ) {
          // This block is not adjacent to the current node
          if ( test_and_apply(to, isolated_block_gain, true /* no tie breaking */ ) ) {
            non_adjacent_block.push_back(to);
          }
        }
      }

      if ( non_adjacent_block.size() > 0 ) {
        // Choose one at random
        const PartitionID to = non_adjacent_block[
          rand.getRandomInt(0, static_cast<int>(non_adjacent_block.size() - 1), cpu_id)];
        best_move.to = to;
        best_move.gain = isolated_block_gain;
      }
    }

    tmp_scores.clear();
    return best_move;
  }


  template<typename PartitionedHypergraph>
  Move_with_transformed_gain computeMaxGainMove_with_transformed_gains(const PartitionedHypergraph& phg,
                          const HypernodeID hn,
                          const bool rebalance = false,
                          const bool consider_non_adjacent_blocks = false,
                          const bool allow_imbalance = false) {
    Derived* derived = static_cast<Derived*>(this);
    RatingMap& tmp_scores = _tmp_scores.local();
    Gain isolated_block_gain = 0;
    derived->precomputeGains(phg, hn, tmp_scores, isolated_block_gain, consider_non_adjacent_blocks);

    auto balance_gain = [&](const PartitionedHypergraph& phg, HypernodeID node, PartitionID from, PartitionID to){
      double gain = 0.0;
      for(int i = 0; i < dimension; i++){
        gain += std::max(0, std::min(phg.nodeWeight(node).weights[i], phg.partWeight(to).weights[i] + 
        phg.nodeWeight(node).weights[i] - _context.partition.max_part_weights[to].weights[i])) * _context.partition.max_part_weights_inv[to][i]
        - std::max(0, std::min(phg.nodeWeight(node).weights[i], phg.partWeight(from).weights[i] - _context.partition.max_part_weights[from].weights[i])) * _context.partition.max_part_weights_inv[from][i];
      }
      return gain;
    };

    PartitionID from = phg.partID(hn);
    Move_with_transformed_gain best_move { from, from, hn, rebalance ? std::numeric_limits<double>::max() : 0.0 };
    HypernodeWeight hn_weight = phg.nodeWeight(hn);
    int cpu_id = THREAD_ID;
    utils::Randomize& rand = utils::Randomize::instance();
    auto test_and_apply = [&](const PartitionID to,
                              const Gain score,
                              const bool no_tie_breaking = false) {
      double balance = balance_gain(phg, hn, from, to);
      if(score > 0 && balance == 0){
        return false;
      }
      double score_gain = score > 0 ? -score / balance : -score * balance;
      bool new_best_gain = balance <= 0 && ((score_gain < best_move.gain) ||
                            (score_gain == best_move.gain &&
                            !_disable_randomization &&
                            (no_tie_breaking || rand.flipCoin(cpu_id))));
      if (new_best_gain) {
        best_move.to = to;
        best_move.gain = score;
        best_move.gain = score_gain;
        return true;
      } else {
        return false;
      }
    };

    for ( const auto& entry : tmp_scores ) {
      const PartitionID to = entry.key;
      if (from != to) {
        const Gain score = derived->gain(entry.value, isolated_block_gain);
        test_and_apply(to, score);
      }
    }

    if ( consider_non_adjacent_blocks && best_move.to == from ) {
      // This is important for our rebalancer as the last fallback strategy
      vec<PartitionID> non_adjacent_block;
      for ( PartitionID to = 0; to < _context.partition.k; ++to ) {
        if ( from != to && !tmp_scores.contains(to) ) {
          // This block is not adjacent to the current node
          if ( test_and_apply(to, isolated_block_gain, true /* no tie breaking */ ) ) {
            non_adjacent_block.push_back(to);
          }
        }
      }

      if ( non_adjacent_block.size() > 0 ) {
        // Choose one at random
        const PartitionID to = non_adjacent_block[
          rand.getRandomInt(0, static_cast<int>(non_adjacent_block.size() - 1), cpu_id)];
        best_move.to = to;
        best_move.gain = isolated_block_gain;
      }
    }

    tmp_scores.clear();
    return best_move;
  }





  inline void computeDeltaForHyperedge(const SynchronizedEdgeUpdate& sync_update) {
    _deltas.local() += AttributedGains::gain(sync_update);
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
    for (const Gain& delta : _deltas) {
      overall_delta += delta;
    }
    return overall_delta;
  }

  void reset() {
    for (Gain& delta : _deltas) {
      delta = 0;
    }
  }

  void changeNumberOfBlocks(const PartitionID new_k) {
    ASSERT(new_k == _context.partition.k);
    for ( auto& tmp_score : _tmp_scores ) {
      if ( static_cast<size_t>(new_k) > tmp_score.size() ) {
        tmp_score = RatingMap(new_k);
      }
    }
    static_cast<Derived*>(this)->changeNumberOfBlocksImpl(new_k);
  }

  
  

private:
  RatingMap constructLocalTmpScores() const {
    return RatingMap(_context.partition.k);
  }

 protected:
  const Context& _context;
  const bool _disable_randomization;
  DeltaGain _deltas;
  TmpScores _tmp_scores;
};

}  // namespace mt_kahypar

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

#pragma once

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/i_rebalancer.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/utils/cast.h"

namespace mt_kahypar {
template <typename GraphAndGainTypes>
class MDRebalancer final : public IRebalancer {
 private:
  using PartitionedHypergraph = typename GraphAndGainTypes::PartitionedHypergraph;
  using GainCache = typename GraphAndGainTypes::GainCache;
  using GainCalculator = typename GraphAndGainTypes::GainComputation;
  using MoveID=uint32_t;
  using PQID=uint64_t;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

public:

  explicit MDRebalancer(HypernodeID , const Context& context, GainCache& gain_cache) :
    _context(context),
    _gain_cache(gain_cache),
    _current_k(context.partition.k),
    _gain(context) { }

  explicit MDRebalancer(HypernodeID num_nodes, const Context& context, gain_cache_t gain_cache) :
    MDRebalancer(num_nodes, context, GainCachePtr::cast<GainCache>(gain_cache)) {}

  MDRebalancer(const MDRebalancer&) = delete;
  MDRebalancer(MDRebalancer&&) = delete;

  MDRebalancer & operator= (const MDRebalancer &) = delete;
  MDRebalancer & operator= (MDRebalancer &&) = delete;

  void initializeImpl(mt_kahypar_partitioned_hypergraph_t&) final;

  bool refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                  const vec<HypernodeID>&,
                  Metrics& best_metrics,
                  double) override final;

  bool refineAndOutputMovesImpl(mt_kahypar_partitioned_hypergraph_t&,
                                const vec<HypernodeID>&,
                                vec<vec<Move>>&,
                                Metrics&,
                                const double) override final;

  bool refineAndOutputMovesLinearImpl(mt_kahypar_partitioned_hypergraph_t&,
                                      const vec<HypernodeID>&,
                                      vec<Move>&,
                                      Metrics&,
                                      const double) override final;
  
  bool labelPropagation(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                       vec<vec<Move>>* moves_by_part,
                                                       vec<Move>* moves_linear,
                                                       Metrics& best_metrics, 
                                                       std::vector<HypernodeID> nodes, double allowed_imbalance, bool balance=false);

private:
  bool refineInternal(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                      vec<vec<Move>>* moves_by_part,
                      vec<Move>* moves_linear,
                      Metrics& best_metrics);

  void simple_lp(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                       vec<vec<Move>>* moves_by_part,
                                                       vec<Move>* moves_linear,
                                                       Metrics& best_metrics, 
                                                       std::vector<HypernodeID> nodes, double allowed_imbalance, bool balance, Gain& local_attributed_gain);

  
  bool greedyRefiner(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                      vec<vec<Move>>* moves_by_part,
                      vec<Move>* moves_linear,
                      Metrics& best_metrics,
                      std::vector<HypernodeID> nodes, bool balance, Gain& local_attributed_gain);

  double L1_balance_gain(PartitionedHypergraph* hypergraph,
                          const HypernodeID node,
                          const PartitionID to);
                        
  bool metis_tiebreak(mt_kahypar_partitioned_hypergraph_t& hypergraph,HypernodeID hn, PartitionID p1, PartitionID p2);


  void interleaveMoveSequenceWithRebalancingMoves(
                      mt_kahypar_partitioned_hypergraph_t& hypergraph,
                      const vec<HypernodeWeight>& initialPartWeights,
                      const std::vector<HypernodeWeight>& max_part_weights,
                      vec<Move>& refinement_moves,
                      vec<vec<Move>>& rebalancing_moves_by_part,
                      vec<Move>& move_order);

  void insertMovesToBalanceBlock(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                                        const PartitionID block,
                                                                        const std::vector<HypernodeWeight>& max_part_weights,
                                                                        const vec<vec<Move>>& rebalancing_moves_by_part,
                                                                        MoveID& next_move_index,
                                                                        vec<HypernodeWeight>& current_part_weights,
                                                                        vec<MoveID>& current_rebalancing_move_index,
                                                                        vec<Move>& move_order);
                                                                    

  int get_max_dimension(mt_kahypar_partitioned_hypergraph_t& hypergraph,
      HypernodeID hn);


  void resizeDataStructuresForCurrentK() {
    // If the number of blocks changes, we resize data structures
    // (can happen during deep multilevel partitioning)
    if ( _current_k != _context.partition.k ) {
      _current_k = _context.partition.k;
      _gain.changeNumberOfBlocks(_current_k);
    }
  }

  const Context& _context;
  GainCache& _gain_cache;
  PartitionID _current_k;
  GainCalculator _gain;
};

}  // namespace kahypar
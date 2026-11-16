/*******************************************************************************
 * MIT License
 *
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2023 Nikolai Maas <nikolai.maas@kit.edu>
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

#include "include/mtkahypartypes.h"

#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/metrics.h"

namespace mt_kahypar {

class IRebalancer {

 public:
  virtual ~IRebalancer() = default;

  void initialize(mt_kahypar_partitioned_hypergraph_t& hypergraph) {
    initializeImpl(hypergraph);
  }

  bool rebalance(mt_kahypar_partitioned_hypergraph_t& hypergraph,
              const parallel::scalable_vector<HypernodeID>& rebalancement_nodes,
              Metrics& best_metrics,
              const bool full_effort) {
    return rebalanceImpl(hypergraph, rebalancement_nodes, best_metrics, full_effort);
  }

  bool rebalanceAndOutputMoves(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                            const parallel::scalable_vector<HypernodeID>& rebalancement_nodes,
                            parallel::scalable_vector<parallel::scalable_vector<Move>>& moves_by_part,
                            Metrics& best_metrics,
                            const bool full_effort) {
    return rebalanceAndOutputMovesImpl(hypergraph, rebalancement_nodes, moves_by_part, best_metrics, full_effort);
  }

  bool rebalanceAndOutputMovesLinear(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                  const parallel::scalable_vector<HypernodeID>& rebalancement_nodes,
                                  parallel::scalable_vector<Move>& moves,
                                  Metrics& best_metrics,
                                  const bool full_effort) {
    return rebalanceAndOutputMovesLinearImpl(hypergraph, rebalancement_nodes, moves, best_metrics, full_effort);
  }

  virtual const HypernodeWeightArray& targetPartWeights() const = 0;

 protected:
  IRebalancer() = default;

 private:
  virtual void initializeImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph) = 0;

  virtual bool rebalanceImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                          const parallel::scalable_vector<HypernodeID>& rebalancement_nodes,
                          Metrics& best_metrics,
                          const bool full_effort) = 0;

  virtual bool rebalanceAndOutputMovesImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                        const parallel::scalable_vector<HypernodeID>& rebalancement_nodes,
                                        parallel::scalable_vector<parallel::scalable_vector<Move>>& moves_by_part,
                                        Metrics& best_metrics,
                                        const bool full_effort) = 0;

  virtual bool rebalanceAndOutputMovesLinearImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                              const parallel::scalable_vector<HypernodeID>& rebalancement_nodes,
                                              parallel::scalable_vector<Move>& moves,
                                              Metrics& best_metrics,
                                              const bool full_effort) = 0;
};

}  // namespace mt_kahypar

/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
******************************************************************************/

#pragma once

#include "tbb/task.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/initial_partitioning/flat/initial_partitioning_data_container.h"

namespace mt_kahypar {

class LabelPropagationInitialPartitioner : public tbb::task {

  using DeltaFunction = std::function<void (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID)>;
  #define NOOP_FUNC [] (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID) { }

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

 public:
  struct MaxGainMove {
    const PartitionID block;
    const Gain gain;
  };

  LabelPropagationInitialPartitioner(const InitialPartitioningAlgorithm,
                                      InitialPartitioningDataContainer& ip_data,
                                      const Context& context,
                                      const int seed) :
    _ip_data(ip_data),
    _context(context),
    _valid_blocks(context.partition.k),
    _tmp_scores(context.partition.k),
    _rng(seed) { }

  tbb::task* execute() override ;

 private:
  bool fitsIntoBlock(PartitionedHypergraph& hypergraph,
                     const HypernodeID hn,
                     const PartitionID block) const {
    ASSERT(block != kInvalidPartition && block < _context.partition.k);
    return hypergraph.partWeight(block) + hypergraph.nodeWeight(hn) <=
      _context.partition.perfect_balance_part_weights[block] *
      std::min(1.005, 1 + _context.partition.epsilon);
  }

  MaxGainMove computeMaxGainMove(PartitionedHypergraph& hypergraph,
                                 const HypernodeID hn) {
    if ( hypergraph.partID(hn) == kInvalidPartition ) {
      return computeMaxGainMoveForUnassignedVertex(hypergraph, hn);
    } else {
      return computeMaxGainMoveForAssignedVertex(hypergraph, hn);
    }
  }

  MaxGainMove computeMaxGainMoveForUnassignedVertex(PartitionedHypergraph& hypergraph,
                                                    const HypernodeID hn);

  MaxGainMove computeMaxGainMoveForAssignedVertex(PartitionedHypergraph& hypergraph,
                                                  const HypernodeID hn);

  MaxGainMove findMaxGainMove(PartitionedHypergraph& hypergraph,
                              const HypernodeID hn,
                              const HypernodeWeight internal_weight);

  void extendBlockToInitialBlockSize(PartitionedHypergraph& hypergraph,
                                     HypernodeID seed_vertex,
                                     const PartitionID block);

  void assignVertexToBlockWithMinimumWeight(PartitionedHypergraph& hypergraph,
                                            const HypernodeID hn);

  InitialPartitioningDataContainer& _ip_data;
  const Context& _context;
  kahypar::ds::FastResetFlagArray<> _valid_blocks;
  parallel::scalable_vector<Gain> _tmp_scores;
  std::mt19937 _rng;
};


} // namespace mt_kahypar

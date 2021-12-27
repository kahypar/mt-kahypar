/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/factories.h"

#include "mt-kahypar/partition/initial_partitioning/flat/random_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/flat/bfs_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/flat/greedy_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/flat/label_propagation_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/flat/policies/gain_computation_policy.h"
#include "mt-kahypar/partition/initial_partitioning/flat/policies/pq_selection_policy.h"

#include "kahypar/meta/abstract_factory.h"
#include "kahypar/meta/registrar.h"

namespace mt_kahypar {
using FlatInitialPartitionerFactory = kahypar::meta::Factory<InitialPartitioningAlgorithm,
        tbb::task *(*)(tbb::task *, const InitialPartitioningAlgorithm, InitialPartitioningDataContainer&,
                       const Context&, const int, const int)>;
}

#define REGISTER_FLAT_INITIAL_PARTITIONER(id, partitioner)                                                     \
  static kahypar::meta::Registrar<FlatInitialPartitionerFactory> register_ ## partitioner(                     \
    id,                                                                                                        \
    [](tbb::task* parent, const InitialPartitioningAlgorithm algorithm,                                        \
       InitialPartitioningDataContainer& ip_hypergraph, const Context& context, const int seed, const int tag) \
    -> tbb::task* {                                                                                            \
    return new(parent->allocate_child()) partitioner(algorithm, ip_hypergraph, context, seed, tag);            \
  })

namespace mt_kahypar {

using GreedyRoundRobinFMInitialPartitioner = GreedyInitialPartitioner<CutGainPolicy, RoundRobinPQSelectionPolicy>;
using GreedyGlobalFMInitialPartitioner = GreedyInitialPartitioner<CutGainPolicy, GlobalPQSelectionPolicy>;
using GreedySequentialFMInitialPartitioner = GreedyInitialPartitioner<CutGainPolicy, SequentialPQSelectionPolicy>;
using GreedyRoundRobinMaxNetInitialPartitioner = GreedyInitialPartitioner<MaxNetGainPolicy, RoundRobinPQSelectionPolicy>;
using GreedyGlobalMaxNetInitialPartitioner = GreedyInitialPartitioner<MaxNetGainPolicy, GlobalPQSelectionPolicy>;
using GreedySequentialMaxNetInitialPartitioner = GreedyInitialPartitioner<MaxNetGainPolicy, SequentialPQSelectionPolicy>;

REGISTER_FLAT_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::random, RandomInitialPartitioner);
REGISTER_FLAT_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::bfs, BFSInitialPartitioner);
REGISTER_FLAT_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::greedy_round_robin_fm, GreedyRoundRobinFMInitialPartitioner);
REGISTER_FLAT_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::greedy_global_fm, GreedyGlobalFMInitialPartitioner);
REGISTER_FLAT_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::greedy_sequential_fm, GreedySequentialFMInitialPartitioner);
REGISTER_FLAT_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::greedy_round_robin_max_net, GreedyRoundRobinMaxNetInitialPartitioner);
REGISTER_FLAT_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::greedy_global_max_net, GreedyGlobalMaxNetInitialPartitioner);
REGISTER_FLAT_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::greedy_sequential_max_net, GreedySequentialMaxNetInitialPartitioner);
REGISTER_FLAT_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::label_propagation, LabelPropagationInitialPartitioner);
}  // namespace mt_kahypar

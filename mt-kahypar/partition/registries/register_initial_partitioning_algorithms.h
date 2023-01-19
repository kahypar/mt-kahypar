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

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/factories.h"

#include "mt-kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/random_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/bfs_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/greedy_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/label_propagation_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/policies/gain_computation_policy.h"
#include "mt-kahypar/partition/initial_partitioning/policies/pq_selection_policy.h"

#include "kahypar/meta/abstract_factory.h"
#include "kahypar/meta/registrar.h"

namespace mt_kahypar {
using InitialPartitionerFactory = kahypar::meta::Factory<InitialPartitioningAlgorithm,
  IInitialPartitioner* (*)(const InitialPartitioningAlgorithm, InitialPartitioningDataContainer&, const Context&, const int, const int)>;
}

#define REGISTER_INITIAL_PARTITIONER(id, partitioner)                                                          \
  static kahypar::meta::Registrar<InitialPartitionerFactory> register_ ## partitioner(                         \
    id,                                                                                                        \
    [](const InitialPartitioningAlgorithm algorithm, InitialPartitioningDataContainer& ip_data,                \
       const Context& context, const int seed, const int tag)                                                  \
    -> IInitialPartitioner* {                                                                                  \
    return new partitioner(algorithm, ip_data, context, seed, tag);                                      \
  })

namespace mt_kahypar {

using GreedyRoundRobinFMInitialPartitioner = GreedyInitialPartitioner<CutGainPolicy, RoundRobinPQSelectionPolicy>;
using GreedyGlobalFMInitialPartitioner = GreedyInitialPartitioner<CutGainPolicy, GlobalPQSelectionPolicy>;
using GreedySequentialFMInitialPartitioner = GreedyInitialPartitioner<CutGainPolicy, SequentialPQSelectionPolicy>;
using GreedyRoundRobinMaxNetInitialPartitioner = GreedyInitialPartitioner<MaxNetGainPolicy, RoundRobinPQSelectionPolicy>;
using GreedyGlobalMaxNetInitialPartitioner = GreedyInitialPartitioner<MaxNetGainPolicy, GlobalPQSelectionPolicy>;
using GreedySequentialMaxNetInitialPartitioner = GreedyInitialPartitioner<MaxNetGainPolicy, SequentialPQSelectionPolicy>;

REGISTER_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::random, RandomInitialPartitioner);
REGISTER_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::bfs, BFSInitialPartitioner);
REGISTER_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::greedy_round_robin_fm, GreedyRoundRobinFMInitialPartitioner);
REGISTER_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::greedy_global_fm, GreedyGlobalFMInitialPartitioner);
REGISTER_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::greedy_sequential_fm, GreedySequentialFMInitialPartitioner);
REGISTER_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::greedy_round_robin_max_net, GreedyRoundRobinMaxNetInitialPartitioner);
REGISTER_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::greedy_global_max_net, GreedyGlobalMaxNetInitialPartitioner);
REGISTER_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::greedy_sequential_max_net, GreedySequentialMaxNetInitialPartitioner);
REGISTER_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::label_propagation, LabelPropagationInitialPartitioner);
}  // namespace mt_kahypar

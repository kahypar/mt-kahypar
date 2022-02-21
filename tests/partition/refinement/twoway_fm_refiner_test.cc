/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
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

#include "gmock/gmock.h"

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/initial_partitioning/flat/bfs_initial_partitioner.h"
#include "mt-kahypar/partition/refinement/fm/sequential_twoway_fm_refiner.h"
#include "mt-kahypar/utils/randomize.h"

using ::testing::Test;

namespace mt_kahypar {

class ATwoWayFmRefiner : public Test {

 public:
  ATwoWayFmRefiner() :
    hypergraph(),
    partitioned_hypergraph(),
    context(),
    refiner(nullptr),
    metrics() {
    context.partition.mode = Mode::direct;
    context.partition.objective = kahypar::Objective::cut;
    context.partition.epsilon = 0.25;
    context.partition.k = 2;
    context.partition.verbose_output = false;

    // Shared Memory
    context.shared_memory.num_threads = 1;

    // Initial Partitioning
    context.initial_partitioning.mode = Mode::deep_multilevel;
    context.initial_partitioning.runs = 1;

    // Read hypergraph
    hypergraph = io::readHypergraphFile(
      "../tests/instances/contracted_ibm01.hgr");
    partitioned_hypergraph = PartitionedHypergraph(
      context.partition.k, hypergraph, parallel_tag_t());
    context.setupPartWeights(hypergraph.totalWeight());
    initialPartition();

    refiner = std::make_unique<SequentialTwoWayFmRefiner>(partitioned_hypergraph, context);
    prng = std::make_unique<std::mt19937>(420);
  }

  void initialPartition() {
    Context ip_context(context);
    ip_context.refinement.label_propagation.algorithm = LabelPropagationAlgorithm::do_nothing;
    InitialPartitioningDataContainer ip_data(partitioned_hypergraph, ip_context);
    BFSInitialPartitioner& initial_partitioner = *new(tbb::task::allocate_root())
      BFSInitialPartitioner(InitialPartitioningAlgorithm::bfs, ip_data, ip_context, 420, 0);
    tbb::task::spawn_root_and_wait(initial_partitioner);
    ip_data.apply();
    metrics.km1 = metrics::km1(partitioned_hypergraph);
    metrics.cut = metrics::hyperedgeCut(partitioned_hypergraph);
    metrics.imbalance = metrics::imbalance(partitioned_hypergraph, context);
  }

  Hypergraph hypergraph;
  PartitionedHypergraph partitioned_hypergraph;
  Context context;
  std::unique_ptr<SequentialTwoWayFmRefiner> refiner;
  Metrics metrics;
  std::unique_ptr<std::mt19937> prng;
};

static constexpr double EPS = 0.05;

TEST_F(ATwoWayFmRefiner, UpdatesImbalanceCorrectly) {
  refiner->refine(metrics, *prng);
  ASSERT_DOUBLE_EQ(metrics::imbalance(partitioned_hypergraph, context), metrics.imbalance);
}

TEST_F(ATwoWayFmRefiner, DoesNotViolateBalanceConstraint) {
  refiner->refine(metrics, *prng);
  ASSERT_LE(metrics.imbalance, context.partition.epsilon + EPS);
}

TEST_F(ATwoWayFmRefiner, UpdatesMetricsCorrectly) {
  refiner->refine(metrics, *prng);
  ASSERT_EQ(metrics::objective(partitioned_hypergraph, context.partition.objective),
            metrics.getMetric(Mode::direct, context.partition.objective));
}

TEST_F(ATwoWayFmRefiner, DoesNotWorsenSolutionQuality) {
  HyperedgeWeight objective_before = metrics::objective(partitioned_hypergraph, context.partition.objective);
    refiner->refine(metrics, *prng);
  ASSERT_LE(metrics.getMetric(Mode::direct, context.partition.objective), objective_before);
}

TEST_F(ATwoWayFmRefiner, DoesProduceCorrectMetricIfExecutedSeveralTimes) {
  refiner->refine(metrics, *prng);
  ASSERT_EQ(metrics::objective(partitioned_hypergraph, context.partition.objective),
            metrics.getMetric(Mode::direct, context.partition.objective));

  partitioned_hypergraph.resetPartition();
  initialPartition();
  refiner->refine(metrics, *prng);
  ASSERT_EQ(metrics::objective(partitioned_hypergraph, context.partition.objective),
            metrics.getMetric(Mode::direct, context.partition.objective));
}
}  // namespace mt_kahypar

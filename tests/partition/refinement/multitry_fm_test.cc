/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "gmock/gmock.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/io/hypergraph_io.h"

#include "mt-kahypar/partition/refinement/fm/multitry_kway_fm.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_strategy.h"

#include "mt-kahypar/partition/initial_partitioning/bfs_initial_partitioner.h"

using ::testing::Test;

namespace mt_kahypar {

class MultiTryFMTest : public ::testing::TestWithParam<PartitionID> {
  using Refiner = MultiTryKWayFM<GainCacheStrategy>;

  public:
    MultiTryFMTest() :
            hypergraph(),
            partitioned_hypergraph(),
            context(),
            refiner(nullptr),
            metrics() {
      context.partition.graph_filename = "../tests/instances/contracted_ibm01.hgr";
      context.partition.graph_community_filename = "../tests/instances/contracted_ibm01.hgr.community";
      context.partition.mode = Mode::direct;
      context.partition.epsilon = 0.25;
      context.partition.verbose_output = false;

      // Shared Memory
      context.shared_memory.num_threads = std::thread::hardware_concurrency();

      // Initial Partitioning
      context.initial_partitioning.mode = Mode::deep_multilevel;
      context.initial_partitioning.runs = 1;

      context.partition.k = GetParam();

      context.refinement.fm.algorithm = FMAlgorithm::fm_gain_cache;
      context.refinement.fm.multitry_rounds = 10;
      context.refinement.fm.num_seed_nodes = 5;

      context.partition.objective = Objective::km1;

      // Read hypergraph
      hypergraph = io::readHypergraphFile(
              "../tests/instances/contracted_unweighted_ibm01.hgr");
      partitioned_hypergraph = PartitionedHypergraph(
              context.partition.k, hypergraph, parallel_tag_t());
      context.setupPartWeights(hypergraph.totalWeight());
      initialPartition();

      refiner = std::make_unique<Refiner>(hypergraph, context);
      refiner->initialize(partitioned_hypergraph);
    }

    void initialPartition() {
      Context ip_context(context);
      ip_context.refinement.label_propagation.algorithm = LabelPropagationAlgorithm::do_nothing;
      InitialPartitioningDataContainer ip_data(partitioned_hypergraph, ip_context);
      BFSInitialPartitioner initial_partitioner(InitialPartitioningAlgorithm::bfs, ip_data, ip_context, 420, 0);
      initial_partitioner.partition();
      ip_data.apply();
      metrics.km1 = metrics::km1(partitioned_hypergraph);
      metrics.cut = metrics::hyperedgeCut(partitioned_hypergraph);
      metrics.imbalance = metrics::imbalance(partitioned_hypergraph, context);
    }

    Hypergraph hypergraph;
    PartitionedHypergraph partitioned_hypergraph;
    Context context;
    std::unique_ptr<Refiner> refiner;
    Metrics metrics;
  };


  TEST_P(MultiTryFMTest, UpdatesImbalanceCorrectly) {
    this->refiner->refine(this->partitioned_hypergraph, {}, this->metrics, std::numeric_limits<double>::max());
    ASSERT_DOUBLE_EQ(metrics::imbalance(this->partitioned_hypergraph, this->context), this->metrics.imbalance);
  }


  TEST_P(MultiTryFMTest, DoesNotViolateBalanceConstraint) {
    this->refiner->refine(this->partitioned_hypergraph, {}, this->metrics, std::numeric_limits<double>::max());
    ASSERT_LE(this->metrics.imbalance, this->context.partition.epsilon);
  }

  TEST_P(MultiTryFMTest, UpdatesMetricsCorrectly) {
    this->refiner->refine(this->partitioned_hypergraph, {}, this->metrics, std::numeric_limits<double>::max());
    ASSERT_EQ(metrics::objective(this->partitioned_hypergraph, this->context.partition.objective),
              this->metrics.getMetric(Mode::direct, this->context.partition.objective));
  }

  TEST_P(MultiTryFMTest, DoesNotWorsenSolutionQuality) {
    HyperedgeWeight objective_before = metrics::objective(this->partitioned_hypergraph, this->context.partition.objective);
    this->refiner->refine(this->partitioned_hypergraph, {}, this->metrics, std::numeric_limits<double>::max());
    ASSERT_LE(this->metrics.getMetric(Mode::direct, this->context.partition.objective), objective_before);
  }

  TEST_P(MultiTryFMTest, AlsoWorksWithNonDefaultFeatures) {
    context.refinement.fm.obey_minimal_parallelism = true;
    context.refinement.fm.rollback_parallel = false;
    context.refinement.fm.perform_moves_global = true;
    HyperedgeWeight objective_before = metrics::objective(this->partitioned_hypergraph, this->context.partition.objective);
    this->refiner->refine(this->partitioned_hypergraph, {}, this->metrics, std::numeric_limits<double>::max());
    ASSERT_LE(this->metrics.getMetric(Mode::direct, this->context.partition.objective), objective_before);
    ASSERT_EQ(metrics::objective(this->partitioned_hypergraph, this->context.partition.objective),
              this->metrics.getMetric(Mode::direct, this->context.partition.objective));
    ASSERT_LE(this->metrics.imbalance, this->context.partition.epsilon);
    ASSERT_DOUBLE_EQ(metrics::imbalance(this->partitioned_hypergraph, this->context), this->metrics.imbalance);
  }

  TEST_P(MultiTryFMTest, WorksWithRefinementNodes) {
    parallel::scalable_vector<HypernodeID> refinement_nodes;
    for (HypernodeID u = 0; u < this->partitioned_hypergraph.initialNumNodes(); ++u) {
      refinement_nodes.push_back(u);
    }
    HyperedgeWeight objective_before = metrics::objective(this->partitioned_hypergraph, this->context.partition.objective);
    this->refiner->refine(this->partitioned_hypergraph, refinement_nodes, this->metrics, std::numeric_limits<double>::max());
    ASSERT_LE(this->metrics.getMetric(Mode::direct, this->context.partition.objective), objective_before);
    ASSERT_EQ(metrics::objective(this->partitioned_hypergraph, this->context.partition.objective),
              this->metrics.getMetric(Mode::direct, this->context.partition.objective));
    ASSERT_LE(this->metrics.imbalance, this->context.partition.epsilon);
    ASSERT_DOUBLE_EQ(metrics::imbalance(this->partitioned_hypergraph, this->context), this->metrics.imbalance);

    std::stringstream buffer;
    std::streambuf* old = std::cout.rdbuf(buffer.rdbuf());  // redirect std::cout to discard output
    this->refiner->printMemoryConsumption();
    std::cout.rdbuf(old);                                   // and reset again
  }


  INSTANTIATE_TEST_CASE_P(
          MultiTryFMTestSuite,
          MultiTryFMTest,
          ::testing::Values(
                  2, 4, 8, 16
          ));




}  // namespace mt_kahypar

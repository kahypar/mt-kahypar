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

#include "gmock/gmock.h"

#include <atomic>

#include "tbb/parallel_invoke.h"

#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/initial_partitioning/flat/random_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/flat/bfs_initial_partitioner.h"
#include "tests/datastructures/hypergraph_fixtures.h"

using ::testing::Test;

namespace mt_kahypar {

template<template<typename> class InitialPartitionerT,
         PartitionID k, size_t runs>
struct TestConfig {
  template<typename TypeTraits>
  using InitialPartitionerTask = InitialPartitionerT<TypeTraits>;
  static constexpr PartitionID K = k;
  static constexpr size_t RUNS = runs;
};

template<typename TypeTraits,
         class InitialPartitionerTask>
class InitialPartitionerRootTaskT : public tbb::task {
  using HyperGraph = typename TypeTraits::HyperGraph;
  using TBB = typename TypeTraits::TBB;
  using InitialPartitioningDataContainer = InitialPartitioningDataContainerT<TypeTraits>;

 public:
  InitialPartitionerRootTaskT(HyperGraph& hypergraph,
                              const Context& context,
                              const size_t runs) :
    _ip_data(hypergraph, context, TBB::GLOBAL_TASK_GROUP),
    _context(context),
    _runs(runs) { }

  tbb::task* execute() override {
    tbb::task::set_ref_count(_runs + 1);
    for ( size_t i = 0; i < _runs; ++i ) {
      tbb::task::spawn(*new(tbb::task::allocate_child())
        InitialPartitionerTask(_ip_data, _context));
    }
    tbb::task::wait_for_all();
    _ip_data.apply();
    return nullptr;
  }

 private:
  InitialPartitioningDataContainer _ip_data;
  const Context& _context;
  const size_t _runs;
};

template<typename Config>
class AFlatInitialPartitionerTest : public ds::AHypergraph<2> {
 private:
  using Base = ds::AHypergraph<2>;

 public:
  using Base::TypeTraits;
  using Base::HwTopology;
  using Base::TBBArena;
  using Base::TestHypergraph;
  using Base::TestStreamingHypergraph;
  using InitialPartitionerRootTask = InitialPartitionerRootTaskT<
    TypeTraits, typename Config::template InitialPartitionerTask<TypeTraits>>;

  AFlatInitialPartitionerTest() :
    Base(),
    hypergraph(),
    context() {
    context.partition.k = Config::K;
    context.partition.epsilon = 0.2;
    context.partition.objective = kahypar::Objective::km1;
    hypergraph = io::readHypergraphFile<TestHypergraph, TestStreamingHypergraph, TBBArena, HwTopology>(
      "../test_instances/test_instance.hgr", context.partition.k, InitialHyperedgeDistribution::equally);
    for ( const HypernodeID& hn : hypergraph.nodes() ) {
      hypergraph.setCommunityID(hn, 0);
    }
    hypergraph.initializeCommunities();
    context.setupPartWeights(hypergraph.totalWeight());
    utils::Timer::instance().disable();
  }

  void execute() {
    InitialPartitionerRootTask& root_ip_task = *new(tbb::task::allocate_root())
      InitialPartitionerRootTask(hypergraph, context, Config::RUNS);
    tbb::task::spawn_root_and_wait(root_ip_task);
  }

  TestHypergraph hypergraph;
  Context context;
};

typedef ::testing::Types<TestConfig<RandomInitialPartitionerT, 2, 1>,
                         TestConfig<RandomInitialPartitionerT, 2, 2>,
                         TestConfig<RandomInitialPartitionerT, 2, 5>,
                         TestConfig<RandomInitialPartitionerT, 3, 1>,
                         TestConfig<RandomInitialPartitionerT, 3, 2>,
                         TestConfig<RandomInitialPartitionerT, 3, 5>,
                         TestConfig<RandomInitialPartitionerT, 4, 1>,
                         TestConfig<RandomInitialPartitionerT, 4, 2>,
                         TestConfig<RandomInitialPartitionerT, 4, 5>,
                         TestConfig<RandomInitialPartitionerT, 5, 1>,
                         TestConfig<RandomInitialPartitionerT, 5, 2>,
                         TestConfig<RandomInitialPartitionerT, 5, 5>,
                         TestConfig<BFSInitialPartitionerT, 2, 1>,
                         TestConfig<BFSInitialPartitionerT, 2, 2>,
                         TestConfig<BFSInitialPartitionerT, 2, 5>,
                         TestConfig<BFSInitialPartitionerT, 3, 1>,
                         TestConfig<BFSInitialPartitionerT, 3, 2>,
                         TestConfig<BFSInitialPartitionerT, 3, 5>,
                         TestConfig<BFSInitialPartitionerT, 4, 1>,
                         TestConfig<BFSInitialPartitionerT, 4, 2>,
                         TestConfig<BFSInitialPartitionerT, 4, 5>,
                         TestConfig<BFSInitialPartitionerT, 5, 1>,
                         TestConfig<BFSInitialPartitionerT, 5, 2>,
                         TestConfig<BFSInitialPartitionerT, 5, 5> > TestConfigs;

TYPED_TEST_CASE(AFlatInitialPartitionerTest, TestConfigs);

TYPED_TEST(AFlatInitialPartitionerTest, HasValidImbalance) {
  this->execute();

  ASSERT_LE(metrics::imbalance(this->hypergraph, this->context),
            this->context.partition.epsilon);
}

TYPED_TEST(AFlatInitialPartitionerTest, AssginsEachHypernode) {
  this->execute();

  for ( const HypernodeID& hn : this->hypergraph.nodes() ) {
    ASSERT_NE(this->hypergraph.partID(hn), -1);
  }
}

TYPED_TEST(AFlatInitialPartitionerTest, HasNoSignificantLowPartitionWeights) {
  this->execute();

  // Each block should have a weight greater or equal than half the average
  // block weight.
  for ( PartitionID block = 0; block < this->context.partition.k; ++block ) {
    ASSERT_GE(this->hypergraph.partWeight(block),
              this->context.partition.perfect_balance_part_weights[block] / 2);
  }
}


}  // namespace mt_kahypar

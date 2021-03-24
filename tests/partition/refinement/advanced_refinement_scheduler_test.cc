/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/partition/refinement/advanced/advanced_refinement_scheduler.h"

using ::testing::Test;

#define MOVE(HN, FROM, TO) Move { FROM, TO, HN, 0 }

namespace mt_kahypar {

class AAdvancedRefinementScheduler : public Test {
 public:
  AAdvancedRefinementScheduler() :
    hg(HypergraphFactory::construct(TBBNumaArena::GLOBAL_TASK_GROUP,
      7 , 4, { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} }, nullptr, nullptr, true)),
    phg(2, TBBNumaArena::GLOBAL_TASK_GROUP, hg),
    context() {
    context.partition.k = 2;
    context.partition.perfect_balance_part_weights.assign(2, 3);
    context.partition.max_part_weights.assign(2, 4);
    context.partition.objective = kahypar::Objective::km1;

    context.shared_memory.num_threads = 2;
    context.refinement.advanced.algorithm = AdvancedRefinementAlgorithm::mock;
    context.refinement.advanced.num_threads_per_search = 1;
    context.refinement.advanced.num_cut_edges_per_block_pair = 50;
    context.refinement.advanced.max_bfs_distance = 2;

    phg.setOnlyNodePart(0, 0);
    phg.setOnlyNodePart(1, 0);
    phg.setOnlyNodePart(2, 0);
    phg.setOnlyNodePart(3, 0);
    phg.setOnlyNodePart(4, 1);
    phg.setOnlyNodePart(5, 1);
    phg.setOnlyNodePart(6, 1);
    phg.initializePartition(TBBNumaArena::GLOBAL_TASK_GROUP);
  }

  Hypergraph hg;
  PartitionedHypergraph phg;
  Context context;
};

template <class F, class K>
void executeConcurrent(F f1, K f2) {
  std::atomic<int> cnt(0);

  tbb::parallel_invoke([&] {
    cnt++;
    while (cnt < 2) { }
    f1();
  }, [&] {
    cnt++;
    while (cnt < 2) { }
    f2();
  });
}

void verifyPartWeights(const vec<HypernodeWeight> actual_weights,
                       const vec<HypernodeWeight> expected_weights) {
  ASSERT_EQ(actual_weights.size(), expected_weights.size());
  for ( size_t i = 0; i < actual_weights.size(); ++i ) {
    ASSERT_EQ(actual_weights[i], expected_weights[i]);
  }
}

TEST_F(AAdvancedRefinementScheduler, MovesOneVertex) {
  AdvancedRefinementScheduler refiner(hg, context, TBBNumaArena::GLOBAL_TASK_GROUP);
  refiner.initialize(phg);
  MoveSequence sequence { { MOVE(3, 0, 1) }, 1 };

  const HyperedgeWeight improvement = refiner.applyMoves(sequence);
  ASSERT_EQ(sequence.state, MoveSequenceState::SUCCESS);
  ASSERT_EQ(improvement, sequence.expected_improvement);
  ASSERT_EQ(1, phg.partID(3));
  verifyPartWeights(refiner.partWeights(), { 3, 4 });
}

TEST_F(AAdvancedRefinementScheduler, MovesVerticesWithIntermediateBalanceViolation) {
  AdvancedRefinementScheduler refiner(hg, context, TBBNumaArena::GLOBAL_TASK_GROUP);
  refiner.initialize(phg);
  MoveSequence sequence { { MOVE(5, 1, 0), MOVE(1, 0, 1), MOVE(3, 0, 1) }, 1 };

  const HyperedgeWeight improvement = refiner.applyMoves(sequence);
  ASSERT_EQ(sequence.state, MoveSequenceState::SUCCESS);
  ASSERT_EQ(improvement, sequence.expected_improvement);
  ASSERT_EQ(1, phg.partID(1));
  ASSERT_EQ(1, phg.partID(3));
  ASSERT_EQ(0, phg.partID(5));
  verifyPartWeights(refiner.partWeights(), { 3, 4 });
}

TEST_F(AAdvancedRefinementScheduler, MovesAVertexThatWorsenSolutionQuality) {
  AdvancedRefinementScheduler refiner(hg, context, TBBNumaArena::GLOBAL_TASK_GROUP);
  refiner.initialize(phg);
  MoveSequence sequence { { MOVE(0, 0, 1) }, 1 };

  const HyperedgeWeight improvement = refiner.applyMoves(sequence);
  ASSERT_EQ(sequence.state, MoveSequenceState::WORSEN_SOLUTION_QUALITY);
  ASSERT_EQ(improvement, 0);
  ASSERT_EQ(0, phg.partID(0));
  verifyPartWeights(refiner.partWeights(), { 4, 3 });
}

TEST_F(AAdvancedRefinementScheduler, MovesAVertexThatViolatesBalanceConstraint) {
  AdvancedRefinementScheduler refiner(hg, context, TBBNumaArena::GLOBAL_TASK_GROUP);
  refiner.initialize(phg);
  MoveSequence sequence { { MOVE(4, 1, 0) }, 1 };

  const HyperedgeWeight improvement = refiner.applyMoves(sequence);
  ASSERT_EQ(sequence.state, MoveSequenceState::VIOLATES_BALANCE_CONSTRAINT);
  ASSERT_EQ(improvement, 0);
  ASSERT_EQ(1, phg.partID(4));
  verifyPartWeights(refiner.partWeights(), { 4, 3 });
}

TEST_F(AAdvancedRefinementScheduler, MovesTwoVerticesConcurrently) {
  context.partition.max_part_weights.assign(2, 5);
  AdvancedRefinementScheduler refiner(hg, context, TBBNumaArena::GLOBAL_TASK_GROUP);
  refiner.initialize(phg);

  MoveSequence sequence_1 { { MOVE(3, 0, 1) }, 1 };
  MoveSequence sequence_2 { { MOVE(5, 1, 0) }, 0 };
  HypernodeWeight improvement_1 = 0, improvement_2 = 0;
  executeConcurrent([&] {
    improvement_1 = refiner.applyMoves(sequence_1);
    ASSERT_EQ(sequence_1.state, MoveSequenceState::SUCCESS);
    ASSERT_EQ(improvement_1, sequence_1.expected_improvement);
    ASSERT_EQ(1, phg.partID(3));
  }, [&] {
    improvement_2 = refiner.applyMoves(sequence_2);
    ASSERT_EQ(sequence_2.state, MoveSequenceState::SUCCESS);
    ASSERT_EQ(improvement_2, sequence_2.expected_improvement);
    ASSERT_EQ(0, phg.partID(5));
  });

  verifyPartWeights(refiner.partWeights(), { 4, 3 });
}

TEST_F(AAdvancedRefinementScheduler, MovesTwoVerticesConcurrentlyWhereOneViolateBalanceConstraint) {
  AdvancedRefinementScheduler refiner(hg, context, TBBNumaArena::GLOBAL_TASK_GROUP);
  refiner.initialize(phg);

  MoveSequence sequence_1 { { MOVE(3, 0, 1) }, 1 };
  MoveSequence sequence_2 { { MOVE(1, 0, 1) }, 0 };
  HypernodeWeight improvement_1 = 0, improvement_2 = 0;
  executeConcurrent([&] {
    improvement_1 = refiner.applyMoves(sequence_1);
  }, [&] {
    improvement_2 = refiner.applyMoves(sequence_2);
  });

  ASSERT_TRUE(sequence_1.state == MoveSequenceState::VIOLATES_BALANCE_CONSTRAINT ||
              sequence_2.state == MoveSequenceState::VIOLATES_BALANCE_CONSTRAINT);
  ASSERT_TRUE(sequence_1.state == MoveSequenceState::SUCCESS ||
              sequence_2.state == MoveSequenceState::SUCCESS);
  if ( sequence_1.state == MoveSequenceState::SUCCESS ) {
    ASSERT_EQ(improvement_1, sequence_1.expected_improvement);
    ASSERT_EQ(1, phg.partID(3));
    ASSERT_EQ(improvement_2, 0);
    ASSERT_EQ(0, phg.partID(1));
  } else {
    ASSERT_EQ(improvement_1, 0);
    ASSERT_EQ(0, phg.partID(3));
    ASSERT_EQ(improvement_2, sequence_2.expected_improvement);
    ASSERT_EQ(1, phg.partID(1));
  }
  verifyPartWeights(refiner.partWeights(), { 3, 4 });
}
}
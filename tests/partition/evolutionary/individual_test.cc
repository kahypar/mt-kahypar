#include <atomic>
#include <type_traits>
#include <mt-kahypar/parallel/tbb_initializer.h>

#include "gmock/gmock.h"

#include "tests/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/evolutionary/individual.h"


using ::testing::Test;

static_assert(!std::is_copy_constructible<mt_kahypar::Individual>::value,
              "Individual must not be copy constructible");
static_assert(!std::is_copy_assignable<mt_kahypar::Individual>::value,
              "Individual must not be copy assignable");
static_assert(std::is_move_constructible<mt_kahypar::Individual>::value,
              "Individual must be move constructible");
static_assert(std::is_move_assignable<mt_kahypar::Individual>::value,
              "Individual must be move assignable");

namespace mt_kahypar { 

template<typename TypeTraits>
class AIndividual : public Test {

 public:
 using Hypergraph = typename TypeTraits::Hypergraph;
 using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
 using Individual = mt_kahypar::Individual;
 using Factory = typename Hypergraph::Factory;

  AIndividual() :
    hypergraph(Factory::construct(
      7 , 4, { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} })),
    partitioned_hypergraph(3, hypergraph, parallel_tag_t()),
    context(createTestContext()) {
    initializePartition();
}

  void initializePartition() {
    if ( hypergraph.nodeIsEnabled(0) ) partitioned_hypergraph.setNodePart(0, 0);
    if ( hypergraph.nodeIsEnabled(1) ) partitioned_hypergraph.setNodePart(1, 0);
    if ( hypergraph.nodeIsEnabled(2) ) partitioned_hypergraph.setNodePart(2, 0);
    if ( hypergraph.nodeIsEnabled(3) ) partitioned_hypergraph.setNodePart(3, 1);
    if ( hypergraph.nodeIsEnabled(4) ) partitioned_hypergraph.setNodePart(4, 1);
    if ( hypergraph.nodeIsEnabled(5) ) partitioned_hypergraph.setNodePart(5, 2);
    if ( hypergraph.nodeIsEnabled(6) ) partitioned_hypergraph.setNodePart(6, 2);
  }

  static Context createTestContext() {
    Context ctx;
    ctx.partition.objective = Objective::km1;  
    return ctx;
  }

  Individual makeDefault() { return Individual(); }
  Individual makeWithFitness(HyperedgeWeight f) { return Individual(f); }
  Individual makeWithPartition(const std::vector<PartitionID>& p) { return Individual(p); }
  Individual makeFromPHG() { return Individual(partitioned_hypergraph, context); }


  Hypergraph hypergraph;
  PartitionedHypergraph partitioned_hypergraph;
  Context context;

};
TYPED_TEST_SUITE(AIndividual, tests::HypergraphTestTypeTraits);

TYPED_TEST(AIndividual, CanBeDefaultConstructed) {
  auto ind = this->makeDefault();
  SUCCEED();
}

TYPED_TEST(AIndividual, CanBeConstructedWithFitness) {
  HyperedgeWeight fitness = 42;
  auto ind = this->makeWithFitness(fitness);
  ASSERT_EQ(fitness, ind.fitness());
}

TYPED_TEST(AIndividual, CanBeConstructedWithPartition) {
  std::vector<PartitionID> partition = {0, 0, 0, 1, 1, 2, 2};
  auto ind = this->makeWithPartition(partition);
  ASSERT_EQ(partition, ind.partition());
}

TYPED_TEST(AIndividual, DefaultConstructedIndividualDiesOnPartitionAccess) {
  auto ind = this->makeDefault();
  EXPECT_DEBUG_DEATH({ (void) ind.partition(); }, "");
}

TYPED_TEST(AIndividual, DefaultConstructedIndividualDiesOnCutEdgesAccess) {
  auto ind = this->makeDefault();
  EXPECT_DEBUG_DEATH({ (void) ind.cutEdges(); }, "");
}

TYPED_TEST(AIndividual, DefaultConstructedIndividualDiesOnStrongCutEdgesAccess) {
  auto ind = this->makeDefault();
  EXPECT_DEBUG_DEATH({ (void) ind.strongCutEdges(); }, "");
}

TYPED_TEST(AIndividual, PartitionOnlyIndividualDiesOnFitnessAccess) {
  std::vector<PartitionID> partition = {0, 0, 0, 1, 1, 2, 2};
  auto ind = this->makeWithPartition(partition);
  EXPECT_DEBUG_DEATH({ (void) ind.fitness(); }, "");
}

TYPED_TEST(AIndividual, CanBeConstructedFromPartitionedHypergraph) {
  auto ind = this->makeFromPHG();
  
  const auto& partition = ind.partition();
  ASSERT_EQ(partition.size(), 7);
  for (HypernodeID hn : this->hypergraph.nodes()) {
    EXPECT_EQ(partition[hn], this->partitioned_hypergraph.partID(hn));
  }
}

TYPED_TEST(AIndividual, PartitionedHypergraphConstructorComputesFitnessCorrectly) {
  auto ind = this->makeFromPHG();
  HyperedgeWeight expected_fitness = mt_kahypar::metrics::quality(this->partitioned_hypergraph, this->context);
  EXPECT_EQ(ind.fitness(), expected_fitness);
}

TYPED_TEST(AIndividual, PartitionedHypergraphConstructorIdentifiesCutEdgesCorrectly) {
  auto ind = this->makeFromPHG();
  const auto& cut_edges = ind.cutEdges();
  
  // Every entry in cut_edges must have connectivity > 1
  for (HyperedgeID he : cut_edges) {
    EXPECT_GT(this->partitioned_hypergraph.connectivity(he), 1);
  }
  
  // Every edge with connectivity > 1 must be in cut_edges
  for (HyperedgeID he : this->partitioned_hypergraph.edges()) {
    if (this->partitioned_hypergraph.connectivity(he) > 1) {
      auto it = std::find(cut_edges.begin(), cut_edges.end(), he);
      EXPECT_NE(it, cut_edges.end()) 
        << "Edge " << he << " with connectivity " 
        << this->partitioned_hypergraph.connectivity(he) << " missing from cut_edges";
    }
  }
}

TYPED_TEST(AIndividual, PartitionedHypergraphConstructorComputesStrongCutEdgesCorrectly) {
  auto ind = this->makeFromPHG();
  const auto& strong_cut_edges = ind.strongCutEdges();
  
  // Compute expected size: sum of (connectivity - 1) for all cut edges
  size_t expected_size = 0;
  for (HyperedgeID he : this->partitioned_hypergraph.edges()) {
    if (this->partitioned_hypergraph.connectivity(he) > 1) {
      expected_size += (this->partitioned_hypergraph.connectivity(he) - 1);
    }
  }
  
  EXPECT_EQ(strong_cut_edges.size(), expected_size);
  
  // Verify each edge appears (connectivity - 1) times
  for (HyperedgeID he : this->partitioned_hypergraph.edges()) {
    if (this->partitioned_hypergraph.connectivity(he) > 1) {
      size_t count = std::count(strong_cut_edges.begin(), strong_cut_edges.end(), he);
      size_t expected_count = this->partitioned_hypergraph.connectivity(he) - 1;
      EXPECT_EQ(count, expected_count) 
        << "Edge " << he << " appears " << count << " times, expected " << expected_count;
    }
  }
}

TYPED_TEST(AIndividual, CutEdgesAreSorted) {
  auto ind = this->makeFromPHG();
  const auto& cut_edges = ind.cutEdges();
  
  EXPECT_TRUE(std::is_sorted(cut_edges.begin(), cut_edges.end()));
}

TYPED_TEST(AIndividual, StrongCutEdgesAreSorted) {
  auto ind = this->makeFromPHG();
  const auto& strong_cut_edges = ind.strongCutEdges();
  
  EXPECT_TRUE(std::is_sorted(strong_cut_edges.begin(), strong_cut_edges.end()));
}

TYPED_TEST(AIndividual, CopyConstructorCreatesIdenticalIndividual) {
  auto ind = this->makeFromPHG();
  auto copy = ind.copy();
  
  EXPECT_EQ(ind.fitness(), copy.fitness());
  EXPECT_EQ(ind.partition(), copy.partition());
  EXPECT_EQ(ind.cutEdges(), copy.cutEdges());
  EXPECT_EQ(ind.strongCutEdges(), copy.strongCutEdges());
}

TYPED_TEST(AIndividual, MoveConstructorMovesIndividual) {
  auto ind = this->makeFromPHG();
  auto copy = ind.copy();
  
  Individual moved(std::move(ind));
  
  EXPECT_EQ(copy.fitness(), moved.fitness());
  EXPECT_EQ(copy.partition(), moved.partition());
  EXPECT_EQ(copy.cutEdges(), moved.cutEdges());
  EXPECT_EQ(copy.strongCutEdges(), moved.strongCutEdges());
}

TYPED_TEST(AIndividual, MoveAssignmentMovesIndividual) {
  auto ind = this->makeFromPHG();
  auto copy = ind.copy();
  
  Individual moved;
  moved = std::move(ind);
  
  EXPECT_EQ(copy.fitness(), moved.fitness());
  EXPECT_EQ(copy.partition(), moved.partition());
  EXPECT_EQ(copy.cutEdges(), moved.cutEdges());
  EXPECT_EQ(copy.strongCutEdges(), moved.strongCutEdges());
}

TYPED_TEST(AIndividual, PartitionedHypergraphConstructorHandlesHigherConnectivityEdge) {
  using Hypergraph = typename TestFixture::Hypergraph;
  using PartitionedHypergraph = typename TestFixture::PartitionedHypergraph;
  using Factory = typename TestFixture::Factory;

  Hypergraph hypergraph = Factory::construct(4, 1, {{0, 1, 2, 3}});
  PartitionedHypergraph partitioned_hypergraph(4, hypergraph, parallel_tag_t());

  if (hypergraph.nodeIsEnabled(0)) partitioned_hypergraph.setNodePart(0, 0);
  if (hypergraph.nodeIsEnabled(1)) partitioned_hypergraph.setNodePart(1, 1);
  if (hypergraph.nodeIsEnabled(2)) partitioned_hypergraph.setNodePart(2, 2);
  if (hypergraph.nodeIsEnabled(3)) partitioned_hypergraph.setNodePart(3, 3);

  Context context = TestFixture::createTestContext();
  mt_kahypar::Individual ind(partitioned_hypergraph, context);

  ASSERT_EQ(ind.cutEdges().size(), 1);
  EXPECT_EQ(ind.strongCutEdges().size(), 3);
  EXPECT_EQ(std::count(ind.strongCutEdges().begin(), ind.strongCutEdges().end(), 0), 3);
}

TYPED_TEST(AIndividual, PartitionedHypergraphConstructorCanProduceNoCutEdges) {
  using Hypergraph = typename TestFixture::Hypergraph;
  using PartitionedHypergraph = typename TestFixture::PartitionedHypergraph;
  using Factory = typename TestFixture::Factory;

  Hypergraph hypergraph = Factory::construct(4, 2, {{0, 1}, {2, 3}});
  PartitionedHypergraph partitioned_hypergraph(2, hypergraph, parallel_tag_t());

  if (hypergraph.nodeIsEnabled(0)) partitioned_hypergraph.setNodePart(0, 0);
  if (hypergraph.nodeIsEnabled(1)) partitioned_hypergraph.setNodePart(1, 0);
  if (hypergraph.nodeIsEnabled(2)) partitioned_hypergraph.setNodePart(2, 1);
  if (hypergraph.nodeIsEnabled(3)) partitioned_hypergraph.setNodePart(3, 1);

  Context context = TestFixture::createTestContext();
  mt_kahypar::Individual ind(partitioned_hypergraph, context);

  EXPECT_DEBUG_DEATH({ (void) ind.cutEdges(); }, "");
  EXPECT_DEBUG_DEATH({ (void) ind.strongCutEdges(); }, "");
}

TYPED_TEST(AIndividual, PrintMethodsDoNotCrashForInitializedIndividual) {
  auto ind = this->makeFromPHG();
  EXPECT_NO_FATAL_FAILURE(ind.print());
  EXPECT_NO_FATAL_FAILURE(ind.printDebug());
}
}
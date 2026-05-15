#include <atomic>
#include <array>
#include <limits>
#include <random>
#include <set>
#include <sstream>
#include <type_traits>
#include <vector>
#include <thread>
#include <mt-kahypar/parallel/tbb_initializer.h>

#include "gmock/gmock.h"

#include "tests/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/evolutionary/individual.h"
#include "mt-kahypar/partition/evolutionary/population.h"
#include "mt-kahypar/partition/context.h"

using ::testing::Test;

namespace mt_kahypar {

template <class F1, class F2>
void executeConcurrent(const F1& f1, const F2& f2) {
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

template <typename TypeTraits>
class APopulation : public Test {

 public:
    using Population = mt_kahypar::Population;
    using Hypergraph = typename TypeTraits::Hypergraph;
    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
    using Factory = typename Hypergraph::Factory;

    APopulation() :
            context(),
            population() {
        initializeContext();
        initializeFitnessOnlyPopulation();
    }

    void initializeContext() {
        context.evolutionary.population_size = 5;
        context.evolutionary.replace_strategy = EvoReplaceStrategy::worst;
        context.partition.objective = Objective::km1;
    }

    void addStartingFitnessIndividual(Population& pop, Context& ctx, const HyperedgeWeight fitness) {
        pop.addStartingIndividual(std::make_shared<Individual>(fitness), ctx);
    }

    void initializeFitnessOnlyPopulation() {
        addStartingFitnessIndividual(population, context, 2);
        addStartingFitnessIndividual(population, context, 4);
        addStartingFitnessIndividual(population, context, 6);
        addStartingFitnessIndividual(population, context, 8);
        addStartingFitnessIndividual(population, context, 10);
    }

    std::shared_ptr<Individual> makeIndividualFromAssignment(const std::array<PartitionID, 7>& assignment) {
        Hypergraph hg = Factory::construct(
            7, 4, { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} });
        PartitionedHypergraph phg(3, hg, parallel_tag_t());

        for (HypernodeID hn = 0; hn < assignment.size(); ++hn) {
            if (hg.nodeIsEnabled(hn)) {
                phg.setNodePart(hn, assignment[hn]);
            }
        }

        return std::make_shared<Individual>(phg, context);
    }

    void fillStructuredPopulation(Population& pop,
                                  const std::vector<std::array<PartitionID, 7>>& assignments) {
        Context local_context = context;
        local_context.evolutionary.population_size = assignments.size();
        for (const auto& assignment : assignments) {
            auto ind = makeIndividualFromAssignment(assignment);
            pop.addStartingIndividual(ind, local_context);
        }
    }

    bool containsFitness(const Population& pop, const HyperedgeWeight fitness) const {
        for (size_t i = 0; i < pop.size(); ++i) {
            if (pop.individualAt(i)->fitness() == fitness) {
                return true;
            }
        }
        return false;
    }

    bool containsFitnessBy(const Population& pop, const HyperedgeWeight fitness) const {
        // Helper to check if population contains specific fitness
        for (size_t i = 0; i < pop.size(); ++i) {
            if (pop.individualAt(i)->fitness() == fitness) {
                return true;
            }
        }
        return false;
    }

    Context context;
    Population population;
};

TYPED_TEST_SUITE(APopulation, tests::HypergraphTestTypeTraits);

// =============================================================================
// Basic Accessor Tests
// =============================================================================

TYPED_TEST(APopulation, PopulationSizeIsCorrect) {
    EXPECT_EQ(this->population.size(), 5);
}

TYPED_TEST(APopulation, BestFitnessReturnsMinimum) {
    EXPECT_EQ(this->population.bestFitness(), 2);
}

TYPED_TEST(APopulation, BestIndividualHasOptimalFitness) {
    auto best = this->population.bestInd();
    EXPECT_EQ(best->fitness(), 2);
}

TYPED_TEST(APopulation, WorstIndividualHasWorstFitness) {
    auto worst = this->population.worstInd();
    EXPECT_EQ(worst->fitness(), 10);
}

TYPED_TEST(APopulation, BestFitnessOnEmptyPopulationReturnsMaxInt) {
    typename TestFixture::Population pop;
    EXPECT_EQ(pop.bestFitness(), std::numeric_limits<int>::max());
}

TYPED_TEST(APopulation, FitnessAtReturnsFitnessAtPosition) {
    EXPECT_EQ(this->population.fitnessAt(0), 2);
    EXPECT_EQ(this->population.fitnessAt(1), 4);
    EXPECT_EQ(this->population.fitnessAt(4), 10);
}

TYPED_TEST(APopulation, IndividualAtReturnsSharedPtr) {
    auto ind = this->population.individualAt(0);
    EXPECT_NE(ind, nullptr);
    EXPECT_EQ(ind->fitness(), 2);
}

// =============================================================================
// Population Modification Tests
// =============================================================================

TYPED_TEST(APopulation, InsertTooManyStartingIndividualsDies) {
    typename TestFixture::Population pop;
    Context ctx = this->context;
    ctx.evolutionary.population_size = 3;

    pop.addStartingIndividual(std::make_shared<Individual>(1), ctx);
    pop.addStartingIndividual(std::make_shared<Individual>(2), ctx);
    pop.addStartingIndividual(std::make_shared<Individual>(3), ctx);

    // Adding more than population_size should trigger assertion
    EXPECT_DEBUG_DEATH({
        pop.addStartingIndividual(std::make_shared<Individual>(4), ctx);
    }, "");
}

TYPED_TEST(APopulation, InsertsIndividualAtWorstPositionForWorstStrategy) {
    this->context.evolutionary.replace_strategy = EvoReplaceStrategy::worst;
    auto ind_good = std::make_shared<Individual>(1);
    auto ind_bad = std::make_shared<Individual>(12);

    const size_t pos1 = this->population.insert(std::move(ind_bad), this->context);
    EXPECT_EQ(pos1, 4);  // Replaces worst (fitness 10)
    EXPECT_EQ(this->population.size(), 5);

    const size_t pos2 = this->population.insert(std::move(ind_good), this->context);
    EXPECT_EQ(pos2, 4);  // Replaces worst again
    EXPECT_EQ(this->population.size(), 5);
}

TYPED_TEST(APopulation, DiverseReplacementPrefersMostSimilarCandidate) {
    typename TestFixture::Population pop;
    this->fillStructuredPopulation(pop, {
        {0, 0, 0, 1, 1, 2, 2},
        {0, 1, 2, 0, 1, 2, 0},
        {0, 1, 2, 2, 1, 2, 0}
    });
    this->context.evolutionary.replace_strategy = EvoReplaceStrategy::diverse;
    auto candidate = this->makeIndividualFromAssignment({0, 1, 2, 0, 1, 2, 0});

    const size_t replaced = pop.insert(std::move(candidate), this->context);
    EXPECT_EQ(replaced, 1);
}

TYPED_TEST(APopulation, StrongDiverseReplacementPrefersMostSimilarCandidate) {
    typename TestFixture::Population pop;
    this->fillStructuredPopulation(pop, {
        {0, 0, 0, 1, 1, 2, 2},
        {0, 1, 2, 0, 1, 2, 0},
        {0, 1, 2, 2, 1, 2, 0}
    });
    this->context.evolutionary.replace_strategy = EvoReplaceStrategy::strong_diverse;
    auto candidate = this->makeIndividualFromAssignment({0, 1, 2, 0, 1, 2, 0});

    const size_t replaced = pop.insert(std::move(candidate), this->context);
    EXPECT_EQ(replaced, 1);
}

TYPED_TEST(APopulation, StrongDiverseRejectsCandidateWorseThanWorstFitness) {
    typename TestFixture::Population pop;
    this->fillStructuredPopulation(pop, {
        {0, 0, 0, 1, 1, 2, 2},
        {0, 1, 2, 0, 1, 2, 0},
        {0, 1, 2, 2, 1, 2, 0}
    });
    this->context.evolutionary.replace_strategy = EvoReplaceStrategy::strong_diverse;
    const HyperedgeWeight old_worst = pop.worstInd()->fitness();
    auto bad_candidate = std::make_shared<Individual>(old_worst + 100);

    const size_t replaced = pop.insert(std::move(bad_candidate), this->context);
    EXPECT_EQ(replaced, std::numeric_limits<unsigned>::max());
    EXPECT_EQ(pop.worstInd()->fitness(), old_worst);
}

TYPED_TEST(APopulation, DiverseRejectsCandidateWorseThanWorstFitness) {
    typename TestFixture::Population pop;
    this->fillStructuredPopulation(pop, {
        {0, 0, 0, 1, 1, 2, 2},
        {0, 1, 2, 0, 1, 2, 0},
        {0, 1, 2, 2, 1, 2, 0}
    });
    this->context.evolutionary.replace_strategy = EvoReplaceStrategy::diverse;
    const HyperedgeWeight old_worst = pop.worstInd()->fitness();
    auto bad_candidate = std::make_shared<Individual>(old_worst + 100);

    const size_t replaced = pop.insert(std::move(bad_candidate), this->context);
    EXPECT_EQ(replaced, std::numeric_limits<unsigned>::max());
    EXPECT_EQ(pop.worstInd()->fitness(), old_worst);
}

// =============================================================================
// List and Selection Tests
// =============================================================================

TYPED_TEST(APopulation, ListOfBestReturnsSortedPrefix) {
    Individuals best3 = this->population.listOfBest(3);
    ASSERT_EQ(best3.size(), 3);
    EXPECT_EQ(best3[0].get()->fitness(), 2);
    EXPECT_EQ(best3[1].get()->fitness(), 4);
    EXPECT_EQ(best3[2].get()->fitness(), 6);
}

TYPED_TEST(APopulation, SampleKParentsReturnBestReturnsKIndividuals) {
    std::mt19937 rng(42);
    std::vector<size_t> parents;
    const size_t k = 3;

    auto best = this->population.sampleKParentsReturnBest(parents, k, true, &rng);

    ASSERT_EQ(parents.size(), k);
    EXPECT_NE(best, nullptr);
    ASSERT_GE(best->fitness(), this->population.bestFitness());
}

TYPED_TEST(APopulation, SampleKParentsReturnsBestAsBestFitness) {
    std::mt19937 rng(42);
    std::vector<size_t> parents;
    const size_t k = 5;

    auto best = this->population.sampleKParentsReturnBest(parents, k, true, &rng);

    EXPECT_EQ(best->fitness(), this->population.bestFitness());
}

TYPED_TEST(APopulation, RandomIndividualReturnsValidIndividual) {
    for (int i = 0; i < 20; ++i) {
        auto ind = this->population.randomIndividual(false);
        EXPECT_NE(ind, nullptr);
        EXPECT_GE(ind->fitness(), this->population.bestFitness());
    }
}

TYPED_TEST(APopulation, DeterministicRandomIndividualIsStable) {
    std::mt19937 rng_a(1337);
    std::mt19937 rng_b(1337);

    auto ind_a = this->population.randomIndividual(true, &rng_a);
    auto ind_b = this->population.randomIndividual(true, &rng_b);

    EXPECT_EQ(ind_a->fitness(), ind_b->fitness());
}
// =============================================================================
// Utility Tests
// =============================================================================

TYPED_TEST(APopulation, ToStringFormatsCsvCorrectly) {
    EXPECT_EQ(Population::toString({}), "");
    EXPECT_EQ(Population::toString({1}), "1");
    EXPECT_EQ(Population::toString({1, 2, 7}), "1,2,7");
    EXPECT_EQ(Population::toString({5, 10, 15, 20}), "5,10,15,20");
}

TYPED_TEST(APopulation, UpdateDiffMatrixReturnsSquareMatrix) {
    typename TestFixture::Population pop;
    this->fillStructuredPopulation(pop, {
        {0, 0, 0, 1, 1, 2, 2},
        {0, 1, 2, 0, 1, 2, 0},
        {0, 1, 2, 2, 1, 2, 0}
    });

    DiffMatrix matrix = pop.updateDiffMatrix();

    ASSERT_EQ(matrix.size(), pop.size());
    for (size_t i = 0; i < matrix.size(); ++i) {
        EXPECT_EQ(matrix[i].size(), pop.size());
        EXPECT_EQ(matrix[i][i], 0);  // Diagonal should be 0
    }
}

TYPED_TEST(APopulation, UpdateDiffMatrixIsSymmetric) {
    typename TestFixture::Population pop;
    this->fillStructuredPopulation(pop, {
        {0, 0, 0, 1, 1, 2, 2},
        {0, 1, 2, 0, 1, 2, 0},
        {0, 1, 2, 2, 1, 2, 0}
    });

    DiffMatrix matrix = pop.updateDiffMatrix();

    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix.size(); ++j) {
            EXPECT_EQ(matrix[i][j], matrix[j][i]) << "Matrix not symmetric at [" << i << ", " << j << "]";
        }
    }
}

TYPED_TEST(APopulation, PrintMethodsDoNotCrash) {
    EXPECT_NO_FATAL_FAILURE(this->population.print());
    EXPECT_NO_FATAL_FAILURE(this->population.printDebug());
}

// =============================================================================
// Thread Safety Tests
// =============================================================================

TYPED_TEST(APopulation, ConcurrentReadsPreserveConsistency) {
    typename TestFixture::Population pop;
    this->fillStructuredPopulation(pop, {
        {0, 0, 0, 1, 1, 2, 2},
        {0, 1, 2, 0, 1, 2, 0},
        {0, 1, 2, 2, 1, 2, 0},
        {0, 1, 0, 1, 1, 2, 2},
        {0, 2, 0, 1, 2, 1, 2}
    });

    const HyperedgeWeight expected_best = pop.bestFitness();

    executeConcurrent([&] {
        for (size_t i = 0; i < 100; ++i) {
            auto ind = pop.randomIndividual(false);
            EXPECT_NE(ind, nullptr);
            EXPECT_EQ(pop.bestFitness(), expected_best);
        }
    }, [&] {
        for (size_t i = 0; i < 100; ++i) {
            auto best = pop.bestInd();
            EXPECT_NE(best, nullptr);
            EXPECT_EQ(best->fitness(), expected_best);
        }
    });
}

TYPED_TEST(APopulation, ConcurrentInsertsPreservePopulationSize) {
    this->context.evolutionary.replace_strategy = EvoReplaceStrategy::worst;

    executeConcurrent([&] {
        for (HyperedgeWeight fitness = 1; fitness <= 200; ++fitness) {
            auto ind = std::make_shared<Individual>(fitness);
            this->population.insert(std::move(ind), this->context);
        }
    }, [&] {
        for (HyperedgeWeight fitness = 201; fitness <= 400; ++fitness) {
            auto ind = std::make_shared<Individual>(fitness);
            this->population.insert(std::move(ind), this->context);
        }
    });

    EXPECT_EQ(this->population.size(), 5);
    EXPECT_EQ(this->population.bestFitness(), 1);
}

TYPED_TEST(APopulation, ConcurrentInsertsAndReadsRemainConsistent) {
    this->context.evolutionary.replace_strategy = EvoReplaceStrategy::worst;

    executeConcurrent([&] {
        for (HyperedgeWeight fitness = 1; fitness <= 150; ++fitness) {
            auto ind = std::make_shared<Individual>(fitness);
            this->population.insert(std::move(ind), this->context);
        }
    }, [&] {
        for (size_t i = 0; i < 150; ++i) {
            auto ind = this->population.randomIndividual(false);
            EXPECT_NE(ind, nullptr);

            HyperedgeWeight fitness = this->population.fitnessAt(0);
            HyperedgeWeight best = this->population.bestFitness();
            EXPECT_LE(best, fitness);
        }
    });

    EXPECT_EQ(this->population.size(), 5);
}

TYPED_TEST(APopulation, DeterministicRandomSampling) {
    std::mt19937 rng_a(1337);
    std::mt19937 rng_b(1337);

    std::vector<size_t> parents_a;
    std::vector<size_t> parents_b;

    auto best_a = this->population.sampleKParentsReturnBest(parents_a, 3, true, &rng_a);
    auto best_b = this->population.sampleKParentsReturnBest(parents_b, 3, true, &rng_b);

    EXPECT_EQ(parents_a, parents_b);
    EXPECT_EQ(best_a->fitness(), best_b->fitness());
}

TYPED_TEST(APopulation, AllPublicMethodsAreThreadSafe) {
    typename TestFixture::Population pop;
    this->fillStructuredPopulation(pop, {
        {0, 0, 0, 1, 1, 2, 2},
        {0, 1, 2, 0, 1, 2, 0},
        {0, 1, 2, 2, 1, 2, 0},
        {0, 1, 0, 1, 1, 2, 2},
        {0, 2, 0, 1, 2, 1, 2}
    });

    // All these should be thread-safe and not deadlock
    executeConcurrent([&] {
        for (int i = 0; i < 50; ++i) {
            pop.size();
            pop.bestFitness();
            pop.bestInd();
            pop.worstInd();
            pop.randomIndividual(false);
            auto partition = pop.bestPartitionCopy();
        }
    }, [&] {
        for (int i = 0; i < 50; ++i) {
            pop.bestFitness();
            pop.fitnessAt(0);
            pop.individualAt(0);
            auto partition = pop.partitionCopyAt(0);
            auto best3 = pop.listOfBest(2);
        }
    });
}

} // namespace mt_kahypar

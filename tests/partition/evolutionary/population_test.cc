#include <atomic>
#include <array>
#include <limits>
#include <set>
#include <sstream>
#include <type_traits>
#include <vector>
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
        Individual ind(fitness);
        pop.addStartingIndividual(ind, ctx);
    }

    void initializeFitnessOnlyPopulation() {
        addStartingFitnessIndividual(population, context, 2);
        addStartingFitnessIndividual(population, context, 4);
        addStartingFitnessIndividual(population, context, 6);
        addStartingFitnessIndividual(population, context, 8);
        addStartingFitnessIndividual(population, context, 10);
    }

    Individual makeIndividualFromAssignment(const std::array<PartitionID, 7>& assignment) {
        Hypergraph hg = Factory::construct(
            7, 4, { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} });
        PartitionedHypergraph phg(3, hg, parallel_tag_t());

        for (HypernodeID hn = 0; hn < assignment.size(); ++hn) {
            if (hg.nodeIsEnabled(hn)) {
                phg.setNodePart(hn, assignment[hn]);
            }
        }

        return Individual(phg, context);
    }

    void fillStructuredPopulation(Population& pop,
                                  const std::vector<std::array<PartitionID, 7>>& assignments) {
        Context local_context = context;
        local_context.evolutionary.population_size = assignments.size();
        for (const auto& assignment : assignments) {
            Individual ind = makeIndividualFromAssignment(assignment);
            pop.addStartingIndividual(ind, local_context);
        }
    }

    bool containsFitness(const Population& pop, const HyperedgeWeight fitness) const {
        for (size_t i = 0; i < pop.size(); ++i) {
            if (pop.individualAt(i).fitness() == fitness) {
                return true;
            }
        }
        return false;
    }

    Context context;
    Population population;
};

TYPED_TEST_SUITE(APopulation, tests::HypergraphTestTypeTraits);

TYPED_TEST(APopulation, InsertTooManyStartingIndividualsDies) {
    Individual ind(15);
    EXPECT_DEBUG_DEATH({ this->population.addStartingIndividual(ind, this->context); }, "");
}

TYPED_TEST(APopulation, AddStartingIndividualReturnsInsertedIndividual) {
    typename TestFixture::Population pop;
    Context local_context = this->context;
    local_context.evolutionary.population_size = 1;

    Individual ind(2);
    const Individual& inserted = pop.addStartingIndividual(ind, local_context);

    EXPECT_EQ(inserted.fitness(), 2);
    EXPECT_EQ(pop.size(), 1);
}

TYPED_TEST(APopulation, BestWorstAndBestFitnessAreComputedCorrectly) {
    EXPECT_EQ(this->population.best(), 0);
    EXPECT_EQ(this->population.worst(), 4);
    EXPECT_EQ(this->population.bestFitness(), 2);
}

TYPED_TEST(APopulation, BestFitnessOnEmptyPopulationReturnsMaxInt) {
    typename TestFixture::Population pop;
    EXPECT_EQ(pop.bestFitness(), std::numeric_limits<int>::max());
}

TYPED_TEST(APopulation, InsertsIndividualAtWorstPositionForWorstStrategy) {
    this->context.evolutionary.replace_strategy = EvoReplaceStrategy::worst;
    Individual ind3(3);
    Individual ind4(6);

    const size_t pos1 = this->population.insert(std::move(ind3), this->context);
    const size_t pos2 = this->population.insert(std::move(ind4), this->context);

    EXPECT_EQ(pos1, 4);
    EXPECT_EQ(pos2, 3);
    EXPECT_EQ(this->population.size(), 5);
}

TYPED_TEST(APopulation, ForceInsertReplacesGivenPosition) {
    Individual ind(5);
    const size_t pos = this->population.forceInsert(std::move(ind), 1);

    EXPECT_EQ(pos, 1);
    EXPECT_EQ(this->population.individualAt(1).fitness(), 5);
}

TYPED_TEST(APopulation, ForceInsertSaveBestPreservesBestIfWorseCandidate) {
    ASSERT_EQ(this->population.best(), 0);
    Individual ind(9);

    const size_t pos = this->population.forceInsertSaveBest(std::move(ind), 0);
    EXPECT_EQ(pos, 0);
    EXPECT_EQ(this->population.individualAt(0).fitness(), 2);
}

TYPED_TEST(APopulation, ForceInsertSaveBestReplacesNonBestEvenIfWorse) {
    ASSERT_EQ(this->population.best(), 0);
    Individual ind(12);
    EXPECT_LE(this->population.individualAt(4).fitness(), 12);

    const size_t pos = this->population.forceInsertSaveBest(std::move(ind), 4);
    EXPECT_EQ(pos, 4);
    EXPECT_EQ(this->population.individualAt(4).fitness(), 12);
}

TYPED_TEST(APopulation, ForceInsertSaveBestReplacesBestIfBetter) {
    ASSERT_EQ(this->population.best(), 0);
    Individual ind(1);

    const size_t pos = this->population.forceInsertSaveBest(std::move(ind), 0);
    EXPECT_EQ(pos, 0);
    EXPECT_EQ(this->population.individualAt(0).fitness(), 1);
}

TYPED_TEST(APopulation, RandomIndividualExceptNeverReturnsException) {
    for (size_t exception = 0; exception < this->population.size(); ++exception) {
        for (size_t i = 0; i < 64; ++i) {
            EXPECT_NE(this->population.randomIndividualExcept(exception), exception);
        }
    }
}

TYPED_TEST(APopulation, RandomIndividualReturnsValidIndex) {
    for (size_t i = 0; i < 128; ++i) {
        const size_t idx = this->population.randomIndividual();
        EXPECT_LT(idx, this->population.size());
    }
}

TYPED_TEST(APopulation, SingleTournamentSelectionReturnsExistingIndividual) {
    const Individual& winner = this->population.singleTournamentSelection();
    EXPECT_TRUE(this->containsFitness(this->population, winner.fitness()));
}

TYPED_TEST(APopulation, TournamentSelectReturnsExistingIndividuals) {
    const auto parents = this->population.tournamentSelect();
    const Individual& first = parents.first.get();
    const Individual& second = parents.second.get();

    EXPECT_TRUE(this->containsFitness(this->population, first.fitness()));
    EXPECT_TRUE(this->containsFitness(this->population, second.fitness()));
}

TYPED_TEST(APopulation, ListOfBestReturnsSortedPrefix) {
    Individuals best3 = this->population.listOfBest(3);
    ASSERT_EQ(best3.size(), 3);
    EXPECT_EQ(best3[0].get().fitness(), 2);
    EXPECT_EQ(best3[1].get().fitness(), 4);
    EXPECT_EQ(best3[2].get().fitness(), 6);
}

TYPED_TEST(APopulation, ToStringFormatsCsvAndEmptyString) {
    EXPECT_EQ(this->population.toString({}), "");
    EXPECT_EQ(this->population.toString({1, 2, 7}), "1,2,7");
}

TYPED_TEST(APopulation, DifferenceUsesCutEdgesAndStrongCutEdges) {
    typename TestFixture::Population pop;
    this->fillStructuredPopulation(pop, {
        {0, 0, 0, 1, 1, 2, 2},
        {0, 1, 2, 0, 1, 2, 0}
    });
    Individual ind = this->makeIndividualFromAssignment({0, 1, 2, 2, 1, 2, 0});

    EXPECT_EQ(pop.difference(ind, 0, false), 1);
    EXPECT_EQ(pop.difference(ind, 0, true), 3);
}

TYPED_TEST(APopulation, DiverseReplacementPrefersMostSimilarCandidate) {
    typename TestFixture::Population pop;
    this->fillStructuredPopulation(pop, {
        {0, 0, 0, 1, 1, 2, 2},
        {0, 1, 2, 0, 1, 2, 0},
        {0, 1, 2, 2, 1, 2, 0}
    });
    this->context.evolutionary.replace_strategy = EvoReplaceStrategy::diverse;
    Individual candidate = this->makeIndividualFromAssignment({0, 1, 2, 0, 1, 2, 0});

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
    Individual candidate = this->makeIndividualFromAssignment({0, 1, 2, 0, 1, 2, 0});

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
    const HyperedgeWeight old_worst = pop.individualAt(pop.worst()).fitness();
    Individual bad_candidate(old_worst + 100);

    const size_t replaced = pop.insert(std::move(bad_candidate), this->context);
    EXPECT_EQ(replaced, std::numeric_limits<unsigned>::max());
    EXPECT_EQ(pop.individualAt(pop.worst()).fitness(), old_worst);
}

TYPED_TEST(APopulation, DiverseRejectsCandidateWorseThanWorstFitness) {
    typename TestFixture::Population pop;
    this->fillStructuredPopulation(pop, {
        {0, 0, 0, 1, 1, 2, 2},
        {0, 1, 2, 0, 1, 2, 0},
        {0, 1, 2, 2, 1, 2, 0}
    });
    this->context.evolutionary.replace_strategy = EvoReplaceStrategy::diverse;
    const HyperedgeWeight old_worst = pop.individualAt(pop.worst()).fitness();
    Individual bad_candidate(old_worst + 100);

    const size_t replaced = pop.insert(std::move(bad_candidate), this->context);
    EXPECT_EQ(replaced, std::numeric_limits<unsigned>::max());
    EXPECT_EQ(pop.individualAt(pop.worst()).fitness(), old_worst);
}

TYPED_TEST(APopulation, UpdateDiffMatrixReturnsSquareCsvWithSeparator) {
    typename TestFixture::Population pop;
    this->fillStructuredPopulation(pop, {
        {0, 0, 0, 1, 1, 2, 2},
        {0, 1, 2, 0, 1, 2, 0},
        {0, 1, 2, 2, 1, 2, 0}
    });
    const std::string matrix = pop.updateDiffMatrix();

    std::stringstream ss(matrix);
    std::string line;
    std::vector<std::string> lines;
    while (std::getline(ss, line)) {
        lines.push_back(line);
    }

    ASSERT_EQ(lines.size(), pop.size() + 1);
    EXPECT_EQ(lines.back(), "---");
    for (size_t i = 0; i < pop.size(); ++i) {
        EXPECT_EQ(std::count(lines[i].begin(), lines[i].end(), ','), pop.size() - 1);
    }
}

TYPED_TEST(APopulation, ConcurrentSafeAccessorsKeepIndicesAndFitnessValid) {
    typename TestFixture::Population pop;
    this->fillStructuredPopulation(pop, {
        {0, 0, 0, 1, 1, 2, 2},
        {0, 1, 2, 0, 1, 2, 0},
        {0, 1, 2, 2, 1, 2, 0},
        {0, 1, 0, 1, 1, 2, 2},
        {0, 2, 0, 1, 2, 1, 2}
    });
    const HyperedgeWeight expected_best_fitness = pop.bestFitnessSafe();

    executeConcurrent([&] {
        for (size_t i = 0; i < 200; ++i) {
            const size_t idx = pop.randomIndividualSafe(this->context);
            EXPECT_LT(idx, pop.size());
            EXPECT_EQ(pop.bestFitnessSafe(), expected_best_fitness);
            EXPECT_LT(pop.bestSafe(), pop.size());
        }
    }, [&] {
        for (size_t i = 0; i < 200; ++i) {
            const size_t idx = pop.randomIndividualSafe(this->context);
            EXPECT_LT(idx, pop.size());
            EXPECT_EQ(pop.fitnessAtSafe(pop.bestSafe()), expected_best_fitness);
            EXPECT_EQ(pop.bestPartitionCopySafe().size(), 7);
            EXPECT_EQ(pop.partitionCopySafe(idx).size(), 7);
            EXPECT_EQ(pop.randomIndividualPartitionCopySafe(this->context).size(), 7);
        }
    });
}

TYPED_TEST(APopulation, IndividualAtSafeMatchesUnsafeAccessor) {
    typename TestFixture::Population pop;
    this->fillStructuredPopulation(pop, {
        {0, 0, 0, 1, 1, 2, 2},
        {0, 1, 2, 0, 1, 2, 0},
        {0, 1, 2, 2, 1, 2, 0},
        {0, 1, 0, 1, 1, 2, 2},
        {0, 2, 0, 1, 2, 1, 2}
    });

    for (size_t i = 0; i < pop.size(); ++i) {
        EXPECT_EQ(pop.individualAtSafe(i).fitness(), pop.individualAt(i).fitness());
    }
}

TYPED_TEST(APopulation, PrintMethodsDoNotCrash) {
    typename TestFixture::Population pop;
    this->fillStructuredPopulation(pop, {
        {0, 0, 0, 1, 1, 2, 2},
        {0, 1, 2, 0, 1, 2, 0},
        {0, 1, 2, 2, 1, 2, 0}
    });

    EXPECT_NO_FATAL_FAILURE(pop.print());
    EXPECT_NO_FATAL_FAILURE(pop.printDebug());
}

TYPED_TEST(APopulation, DeterministicSafeRandomMethodsAreStableAcrossThreads) {
    typename TestFixture::Population pop;
    this->fillStructuredPopulation(pop, {
        {0, 0, 0, 1, 1, 2, 2},
        {0, 1, 2, 0, 1, 2, 0},
        {0, 1, 2, 2, 1, 2, 0},
        {0, 1, 0, 1, 1, 2, 2},
        {0, 2, 0, 1, 2, 1, 2}
    });
    Context ctx(this->context);
    ctx.partition.deterministic = true;
    ctx.partition.seed = 1337;
    std::mt19937 rng(ctx.partition.seed);
    std::vector<size_t> expected_idxs;
    std::vector<std::vector<PartitionID>> expected_partitions;
    expected_partitions.reserve(100);
    for (size_t i = 0; i < 100; ++i) {
        expected_idxs.push_back(pop.randomIndividualSafe(ctx, &rng));
        expected_partitions.push_back(pop.randomIndividualPartitionCopySafe(ctx, &rng));
    }
    executeConcurrent([&] {
        std::mt19937 thread_rng(ctx.partition.seed);
        for (size_t i = 0; i < 100; ++i) {
            EXPECT_EQ(pop.randomIndividualSafe(ctx, &thread_rng), expected_idxs[i]);
            EXPECT_EQ(pop.randomIndividualPartitionCopySafe(ctx, &thread_rng), expected_partitions[i]);
        }
    }, [&] {
        std::mt19937 thread_rng(ctx.partition.seed);
        for (size_t i = 0; i < 100; ++i) {
            EXPECT_EQ(pop.randomIndividualSafe(ctx, &thread_rng), expected_idxs[i]);
            EXPECT_EQ(pop.randomIndividualPartitionCopySafe(ctx, &thread_rng), expected_partitions[i]);
        }
    });
}

TYPED_TEST(APopulation, ConcurrentInsertsPreservePopulationSizeAndBestFitness) {
    this->context.evolutionary.replace_strategy = EvoReplaceStrategy::worst;

    executeConcurrent([&] {
        for (HyperedgeWeight fitness = 1; fitness <= 200; ++fitness) {
            Individual ind(fitness);
            this->population.insert(std::move(ind), this->context);
        }
    }, [&] {
        for (HyperedgeWeight fitness = 201; fitness <= 400; ++fitness) {
            Individual ind(fitness);
            this->population.insert(std::move(ind), this->context);
        }
    });

    EXPECT_EQ(this->population.size(), 5);
    EXPECT_EQ(this->population.bestFitness(), 1);
    EXPECT_LT(this->population.best(), this->population.size());
    EXPECT_LT(this->population.worst(), this->population.size());
}

TYPED_TEST(APopulation, ConcurrentInsertsAndSafeReadsRemainConsistent) {
    this->context.evolutionary.replace_strategy = EvoReplaceStrategy::worst;

    executeConcurrent([&] {
        for (HyperedgeWeight fitness = 1; fitness <= 300; ++fitness) {
            Individual ind(fitness);
            this->population.insert(std::move(ind), this->context);
        }
    }, [&] {
        for (size_t i = 0; i < 300; ++i) {
            const size_t idx = this->population.randomIndividualSafe(this->context);
            EXPECT_LT(idx, this->population.size());
            EXPECT_LT(this->population.bestSafe(), this->population.size());

            const HyperedgeWeight fitness = this->population.fitnessAtSafe(idx);
            const HyperedgeWeight best_fitness = this->population.bestFitnessSafe();

            EXPECT_LE(best_fitness, fitness);
            EXPECT_GE(fitness, 1);
            EXPECT_LE(fitness, 300);
        }
    });

    EXPECT_EQ(this->population.size(), 5);
    EXPECT_LT(this->population.best(), this->population.size());
    EXPECT_LT(this->population.worst(), this->population.size());
}

} // namespace mt_kahypar

/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2017 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2017 Robin Andre <robinandre1995@web.de>
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

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>
#include <mutex>
#include <atomic>

#include "mt-kahypar/partition/evolutionary/individual.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/partitioner.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {
    struct PopulationSummary {
        HyperedgeWeight best_fitness = 0;
        HyperedgeWeight median_fitness = 0;
        double mean_fitness = 0.0;
        HyperedgeWeight worst_fitness = 0;
        size_t min_distance = 0;
        double median_distance = 0.0;
        double mean_distance = 0.0;
        size_t max_distance = 0;
    };

    struct InsertionDetails {
        size_t position = std::numeric_limits<size_t>::max();
        bool was_accepted = false;
        int64_t dist_to_replaced = -1;
        size_t min_dist_pop = 0;
    };

    class Population {
    private:
        static constexpr bool debug = true;

    public:
        explicit Population() : _population_mutex(),
                                _individuals(),
                                _slot_replacements(),
                                total_insert_attempts(0),
                                rejected_worse_than_worst(0),
                                accepted_replacement(0) {
        }
        
        std::atomic<size_t> total_insert_attempts;
        std::atomic<size_t> rejected_worse_than_worst;
        std::atomic<size_t> accepted_replacement;

        // NM: remove distinction between thread-safe / non thread-safe methods, all public methods should be thread-safe

        //Stats accessors
        size_t size() const;
        HyperedgeWeight bestFitness() const;
        //only used to first sort the individuals by fittness and then retrieve the first n individuals (meta evo)
        HyperedgeWeight fitnessAt(size_t pos) const;

        //modifiers
        size_t insert(std::shared_ptr<Individual> individual, const Context &context);
        InsertionDetails insertWithDetails(std::shared_ptr<Individual> individual, const Context &context);
        void addStartingIndividual(std::shared_ptr<Individual> individual, const Context &context);

        //_individuals accessors
        Individuals listOfBest(const size_t &amount) const;
        std::shared_ptr<const Individual> bestInd() const;
        std::shared_ptr<const Individual> worstInd() const;
        std::shared_ptr<const Individual> individualAt(size_t pos) const;
        std::shared_ptr<const Individual> sampleKParentsReturnBest(std::vector<size_t> &parents, size_t k, bool deterministic,
                                                             std::mt19937 *rng = nullptr) const;
        std::shared_ptr<const Individual> randomIndividual(bool deterministic, std::mt19937 *rng = nullptr) const;

        //parition accessors
        std::vector<PartitionID> bestPartitionCopy() const;
        std::vector<PartitionID> randomIndividualPartitionCopy(bool deterministic, std::mt19937 *rng = nullptr);
        std::vector<PartitionID> partitionCopyAt(size_t pos) const;
        std::vector<HyperedgeID> cutEdgesCopyAt(size_t pos) const;

        //summary & diversity
        PopulationSummary computeSummary() const;
        std::vector<size_t> slotReplacements() const;
        static size_t distanceBetween(const Individual& a, const Individual& b, bool strong_set = false);

        //debug
        void print() const;
        void printDebug() const;
        static std::string toString(const std::vector<size_t> &values);
        std::vector<std::vector<size_t>> updateDiffMatrix() const;

    private:
        size_t forceInsert(std::shared_ptr<Individual> individual, size_t position);
        size_t difference(std::shared_ptr<Individual> individual, size_t position, bool strong_set) const;
        size_t diversePosition(const Individuals &individuals, std::shared_ptr<Individual> individual, bool strong_set);
        size_t worst(const Individuals &individuals) const;

        mutable std::mutex _population_mutex;
        std::vector<std::shared_ptr<Individual> > _individuals;
        std::vector<size_t> _slot_replacements;

    };

    std::ostream &operator<<(std::ostream &os, const Population &population);
    using DiffMatrix = std::vector<std::vector<size_t>>;
} // namespace mt_kahypar

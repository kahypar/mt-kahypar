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

#include "mt-kahypar/partition/evolutionary/individual.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/partitioner.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {
class Population {
 private:
  static constexpr bool debug = false;

 public:
  explicit Population() :
    _population_mutex(),
    _individuals() { }

  size_t insert(Individual&& individual, const Context& context);
  size_t forceInsert(Individual&& individual, size_t position) ;
  size_t forceInsertSaveBest(Individual&& individual, size_t position) ;
  const Individual & singleTournamentSelection() const ;
  std::pair<std::reference_wrapper<const Individual>,
                   std::reference_wrapper<const Individual> > tournamentSelect() const ;
  const Individual & addStartingIndividual(Individual& individual, Context& context) ;
  size_t size() const;
  size_t randomIndividual() const ;
  size_t randomIndividualExcept(size_t exception) const ;
  size_t best() const ;
  HyperedgeWeight bestFitness() const ;
  size_t worst() ;
  const Individual & individualAt(size_t pos) const ;

  // Thread-safe accessors
  size_t sampleKParentsReturnBestIndex(std::vector<size_t>& parents, size_t k, bool deterministic, std::mt19937* rng = nullptr) ;
  size_t randomIndividualSafe(bool deterministic, std::mt19937* rng = nullptr) ;
  const Individual& individualAtSafe(size_t pos);
  std::vector<PartitionID> bestPartitionCopySafe() ;
  std::vector<PartitionID> randomIndividualPartitionCopySafe(bool deterministic, std::mt19937* rng = nullptr) ;
  std::vector<PartitionID> partitionCopySafe(size_t pos) ;
  std::vector<HyperedgeID> cutEdgesCopySave(size_t pos);
  HyperedgeWeight bestFitnessSafe() ;
  HyperedgeWeight fitnessAtSafe( size_t pos) ;
  size_t bestSafe() ;
  Individuals listOfBest(const size_t& amount) const ;

  void print() const ;
  void printDebug() const ;
  size_t difference(const Individual& individual,  size_t position, bool strong_set) const ;
  std::string toString(const std::vector<size_t>& values) const ;
  std::string updateDiffMatrix();

 private:
   size_t replaceDiverse(Individual&& individual,  bool strong_set) ;

  std::mutex _population_mutex;
  std::vector<Individual> _individuals;
  std::vector<std::vector<size_t>> _diff_matrix;
  std::mutex _diff_mutex;
};
std::ostream& operator<< (std::ostream& os, const Population& population);
}  // namespace mt_kahypar
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
#include "population.h"

namespace mt_kahypar {
std::ostream& operator<< (std::ostream& os, const Population& population) {
  for (size_t i = 0; i < population.size(); ++i) {
    os << population.individualAt(i)->fitness() << " ";
  }
  return os;
}

size_t Population::insert(std::shared_ptr<Individual> individual, const Context& context) {
    // NM: rewrite this method as follows (requires that Individuals already use shared_ptr):
    // 1. copy the vector of individuals, release lock after copy
    // 2. determine replacement on copied vector
    // 3. take lock and check whether individual at determined index is still the same
    // 4. if yes, replace individual, if no, go back to 1.

    std::lock_guard<std::mutex> guard(_population_mutex);
    DBG << context.evolutionary.replace_strategy;
    switch (context.evolutionary.replace_strategy) {
      case EvoReplaceStrategy::worst:
        return forceInsert(std::move(individual), worst());
      case EvoReplaceStrategy::diverse:
        return replaceDiverse(std::move(individual), false);

      case EvoReplaceStrategy::strong_diverse:
        return replaceDiverse(std::move(individual), true);
      default:
        return std::numeric_limits<int>::max();
    }
  }

  void Population::addStartingIndividual(std::shared_ptr<Individual> individual, const Context& context) {
  std::lock_guard<std::mutex> guard(_population_mutex);
  _individuals.emplace_back(std::move(individual));
    ASSERT(_individuals.size() <= context.evolutionary.population_size);
    DBG << "Individual" << _individuals.size() - 1
        << V(_individuals.back()->fitness());
  }

  size_t Population::forceInsert(std::shared_ptr<Individual> individual, const size_t position) {
    DBG << V(position) << V(individual->fitness());
    _individuals[position] = std::move(individual);

    return position;
  }

  //possibly unnecesary -->
  size_t Population::size() const {
    return _individuals.size();
  }
  size_t Population::randomIndividual() const {
    return utils::Randomize::instance().getRandomInt(0, size() - 1, THREAD_ID);
  }
  size_t Population::randomIndividualExcept(const size_t exception) const {
    size_t target = utils::Randomize::instance().getRandomInt(0, size() - 2, THREAD_ID);
    if (target == exception) {
      target = size() - 1;
    }
    return target;
  }
  // <--

  std::shared_ptr<Individual> Population::best() const {
    std::shared_ptr<Individual> best_individual;
    HyperedgeWeight best_fitness = std::numeric_limits<int>::max();
    std::lock_guard<std::mutex> guard(_population_mutex);
    for (size_t i = 0; i < size(); ++i) {
      const HyperedgeWeight result = _individuals[i]->fitness();
      if (result < best_fitness) {
        best_individual = _individuals[i];
        best_fitness = result;
      }
    }
    DBG << V(best_individual) << V(best_fitness);
    return best_individual;
  }
  std::shared_ptr<Individual> Population::worstInd() const {
  size_t worst_position = std::numeric_limits<size_t>::max();
  HyperedgeWeight worst_fitness = std::numeric_limits<int>::min();
  for (size_t i = 0; i < size(); ++i) {
    HyperedgeWeight result = _individuals[i]->fitness();
    if (result > worst_fitness) {
      worst_position = i;
      worst_fitness = result;
    }
  }
  DBG << V(worst_position) << V(worst_fitness);
  return _individuals[worst_position];
}

  HyperedgeWeight Population::bestFitness() const {
    std::lock_guard<std::mutex> guard(_population_mutex);
    size_t best_position = std::numeric_limits<size_t>::max();
    HyperedgeWeight best_fitness = std::numeric_limits<int>::max();
    if (size() == 0) {
      DBG << "SIZE IS 0";
      return best_fitness;
    }
    for (size_t i = 0; i < size(); ++i) {
      const HyperedgeWeight result = _individuals[i]->fitness();
      if (result < best_fitness) {
        best_position = i;
        best_fitness = result;
      }
    }
    ASSERT(best_position != std::numeric_limits<size_t>::max());
    DBG << V(best_position) << V(best_fitness);
    return best_fitness;
  }


  size_t Population::worst() const {
    size_t worst_position = std::numeric_limits<size_t>::max();
    HyperedgeWeight worst_fitness = std::numeric_limits<int>::min();
    for (size_t i = 0; i < size(); ++i) {
      HyperedgeWeight result = _individuals[i]->fitness();
      if (result > worst_fitness) {
        worst_position = i;
        worst_fitness = result;
      }
    }
    DBG << V(worst_position) << V(worst_fitness);
    return worst_position;
  }

  const std::shared_ptr<Individual> Population::individualAt(const size_t pos) const {
    return _individuals[pos];
  }

  // Thread-safe accessors
  std::shared_ptr<Individual> Population::sampleKParentsReturnBestIndex(std::vector<size_t>& parents, const size_t k, const bool deterministic, std::mt19937* rng) const {
    ASSERT(k > 0);
    ASSERT(k <= _individuals.size());
    std::lock_guard<std::mutex> guard(_population_mutex);
    std::shared_ptr<Individual> best;
    HyperedgeWeight best_fitness = std::numeric_limits<HyperedgeWeight>::max();
    std::vector<size_t> indices(_individuals.size());
    std::iota(indices.begin(), indices.end(), 0);
    //TODO: Make forloop less confusing
    for (size_t t = 0; t < k; ++t) {
      const size_t i = _individuals.size() - 1 - t;
      size_t rnd;
      if (deterministic) {
        ASSERT(rng != nullptr, "Deterministic mode requires a valid RNG");
        std::uniform_int_distribution<size_t> gen(0, i);
        rnd = gen(*rng);
      } else {
        rnd = utils::Randomize::instance().getRandomInt(0, i, THREAD_ID);
      }
      std::swap(indices[i], indices[rnd]);
      auto candidate = indices[i];
      if (const HyperedgeWeight candidate_fitness = _individuals[candidate]->fitness(); candidate_fitness < best_fitness) {
        best_fitness = candidate_fitness;
        best = _individuals[candidate];
      }
      parents.push_back(indices[i]);
    }
    return best;
  }


  std::shared_ptr<Individual> Population::randomIndividualSafe(const bool deterministic, std::mt19937* rng) const {
    std::lock_guard<std::mutex> guard(_population_mutex);
    size_t pos;
    if (deterministic) {
      ASSERT(rng != nullptr, "Deterministic mode requires a valid RNG");
      std::uniform_int_distribution<size_t> dist(0, _individuals.size() - 1);
      pos = dist(*rng);
    } else {
      pos = utils::Randomize::instance().getRandomInt(0, _individuals.size() - 1, THREAD_ID);
    }
    return _individuals[pos];
  }

  const std::shared_ptr<Individual> Population::individualAtSafe(const size_t pos) {
    std::lock_guard<std::mutex> guard(_population_mutex);
    return _individuals[pos];
  }

  //Not thread safe anymore
  std::vector<PartitionID> Population::bestPartitionCopySafe() const {
    return best()->partition();// returns copy
  }

  std::vector<PartitionID> Population::randomIndividualPartitionCopySafe(const bool deterministic, std::mt19937* rng) {
    std::lock_guard<std::mutex> guard(_population_mutex);
    size_t random_position;
    if (deterministic) {
      ASSERT(rng != nullptr, "Deterministic mode requires a valid RNG");
      std::uniform_int_distribution<size_t> dist(0, _individuals.size() - 1);
      random_position = dist(*rng);
    } else {
      random_position = utils::Randomize::instance().getRandomInt(0, _individuals.size() - 1, THREAD_ID);
    }
    return _individuals[random_position]->partition(); // returns copy
  }

  std::vector<PartitionID> Population::partitionCopySafe(const size_t pos) const {
    std::lock_guard<std::mutex> guard(_population_mutex);
    return _individuals[pos]->partition(); // returns copy
  }

  std::vector<HyperedgeID> Population::cutEdgesCopySave(size_t pos) const {
    std::lock_guard<std::mutex> guard(_population_mutex);
    return _individuals[pos]->cutEdges();
  }

  HyperedgeWeight Population::fitnessAtSafe(const size_t pos) const {
    std::lock_guard<std::mutex> guard(_population_mutex);
    return _individuals[pos]->fitness();
  }

  Individuals Population::listOfBest(const size_t& amount) const {
    std::vector<std::pair<HyperedgeWeight, size_t> > sorting;
    for (size_t i = 0; i < _individuals.size(); ++i) {
      sorting.push_back(std::make_pair(_individuals[i]->fitness(), i));
    }

    std::partial_sort(sorting.begin(), sorting.begin() + amount, sorting.end());

    Individuals best_individuals;
    for (size_t i = 0; i < amount; ++i) {
      best_individuals.push_back(_individuals[sorting[i].second]);
    }
    return best_individuals;
  }

  void Population::print() const {
    std::cout << std::endl << "Population Fitness: ";
    for (size_t i = 0; i < _individuals.size(); ++i) {
      std::cout << _individuals[i]->fitness() << " ";
    }
    std::cout << std::endl;
  }
  void Population::printDebug() const {
    for (size_t i = 0; i < _individuals.size(); ++i) {
      _individuals[i]->printDebug();
    }
  }
  size_t Population::difference(std::shared_ptr<Individual> individual, const size_t position,
                           const bool strong_set) const {
    std::vector<HyperedgeID> output_diff;
    if (strong_set) {
      ASSERT(std::is_sorted(_individuals[position]->strongCutEdges().begin(),
                            _individuals[position]->strongCutEdges().end()));
      ASSERT(std::is_sorted(individual->strongCutEdges().begin(),
                            individual->strongCutEdges().end()));
      std::set_symmetric_difference(_individuals[position]->strongCutEdges().begin(),
                                    _individuals[position]->strongCutEdges().end(),
                                    individual->strongCutEdges().begin(),
                                    individual->strongCutEdges().end(),
                                    std::back_inserter(output_diff));
    } else {
      ASSERT(std::is_sorted(_individuals[position]->cutEdges().begin(),
                            _individuals[position]->cutEdges().end()));
      ASSERT(std::is_sorted(individual->cutEdges().begin(),
                            individual->cutEdges().end()));
      std::set_symmetric_difference(_individuals[position]->cutEdges().begin(),
                                    _individuals[position]->cutEdges().end(),
                                    individual->cutEdges().begin(),
                                    individual->cutEdges().end(),
                                    std::back_inserter(output_diff));
    }
    DBG << V(output_diff.size());
    return output_diff.size();
  }

  std::string Population::toString(const std::vector<size_t>& values) const {
    if (values.empty()) {
        return "";
    }

    std::string result;
    for (size_t i = 0; i < values.size(); ++i) {
        result += std::to_string(values[i]);
        if (i < values.size() - 1) {
            result += ",";
        }
    }
    return result;
}

  std::string Population::updateDiffMatrix() {
    std::lock_guard<std::mutex> guard(_diff_mutex);

    // should theoretically only happen once at init
    _diff_matrix.resize(_individuals.size());
    for (size_t i = 0; i < _diff_matrix.size(); ++i) {
      _diff_matrix[i].resize(_individuals.size());
    }

    for (size_t i = 0; i < _individuals.size(); ++i) {
      for (size_t j = 0; j < _individuals.size(); ++j) {
        if (i == j) {
          _diff_matrix[i][j] = 0;
        } else if (i < j) {
          const size_t diff = difference(_individuals[i], j, false);
          _diff_matrix[i][j] = diff;
          _diff_matrix[j][i] = diff;
        }
      }
    }

    std::string matrix_string;
    // return matrix as String
    for (size_t i = 0; i < _diff_matrix.size(); i++) {
      std::string row = toString(_diff_matrix[i]);
      matrix_string += row + "\n";
    }
    // add separator between each matrix
    matrix_string += "---\n";
    return matrix_string;
  }

  size_t Population::replaceDiverse(std::shared_ptr<Individual> individual, const bool strong_set) {
  size_t max_similarity = std::numeric_limits<size_t>::max();
  size_t max_similarity_id = 0;
  if (individual->fitness() > individualAt(worst())->fitness()) {
    DBG << "COLLAPSE";
    return std::numeric_limits<unsigned>::max();
  }
  for (size_t i = 0; i < size(); ++i) {
    if (_individuals[i]->fitness() >= individual->fitness()) {
      const size_t similarity = difference(individual, i, strong_set);
      DBG << "SYMMETRIC DIFFERENCE:" << similarity << " from" << i;
      if (similarity < max_similarity) {
        max_similarity = similarity;
        max_similarity_id = i;
      }
    }
  }
  DBG << V(max_similarity_id) << V(max_similarity);
  forceInsert(std::move(individual), max_similarity_id);
  return max_similarity_id;
  }
}
// namespace mt_kahypar
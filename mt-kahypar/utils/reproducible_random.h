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

#pragma once

#include <random>

#include "mt-kahypar/parallel/parallel_counting_sort.h"
#include "hash.h"

namespace mt_kahypar::utils {

/*!
 * Combines a global seed and an iteration of a loop to initialize an RNG for that iteration
 */
inline size_t seed_iteration(size_t seed, size_t iteration) {
  return hashing::integer::combine(seed, hashing::integer::hash(iteration));
}

template< template<typename> typename UnqualifiedHashFunction >
class UniformRandomSelector {
public:
  using int_type = size_t;
  using HashFunction = UnqualifiedHashFunction<int_type>;
  using RNG = hashing::HashRNG<HashFunction>;

  UniformRandomSelector(const HashFunction& hash_function, int_type seed) : rng(hash_function, seed) { }

  /*!
   * Call when you find an element with the same score as the current best.
   * Selects one of the elements with best score uniformly at random.
   */
  bool replace_sample() {
    return (dist(rng, std::uniform_int_distribution<size_t>::param_type(0, ++counter)) == 0);
  }

  /*!
   * Call when you find an element with better score for the first time
   */
  void replace() {
    counter = 0;
  }

  /*!
   * Does not reseed the hash function but reseeds the RNG, so that it can be used for neighbor rating in coarsening
   */
  void reset(size_t seed) {
    replace();
    rng.init(seed);
  }

private:
  RNG rng;
  std::uniform_int_distribution<int_type> dist;
  size_t counter = 0;
};

class BucketSeeding {
protected:
  void fill_seeds(std::mt19937& rng) {
    for (size_t i = 0; i < num_buckets; ++i) {
      seeds[i] = rng();
    }
  }

  static constexpr size_t num_buckets = 256;
  std::array<std::mt19937::result_type, num_buckets> seeds;
};

template<class T>
class DeterministicParallelUniformRandomShuffler : public BucketSeeding {
public:
  DeterministicParallelUniformRandomShuffler(vec<T>&& external_elements) :
          elements(std::move(external_elements)),
          shuffled_elements(elements.size(), T())
  {

  }

  template<typename F>
  void shuffle(F get_bucket, size_t num_tasks, std::mt19937& rng) {    // TODO implement get_bucket as deterministic pseudorandom number generator
    if (elements.size() < 1 << 17) {
      std::copy_n(elements.begin(), elements.size(), shuffled_elements.begin());
      std::shuffle(shuffled_elements.begin(), shuffled_elements.end(), rng);
      return;
    }

    vec<uint32_t> bucket_bounds = parallel::counting_sort(elements, shuffled_elements, num_buckets, get_bucket, num_tasks);
    ASSERT(bucket_bounds.size() == num_buckets + 1);

    fill_seeds(rng);

    tbb::parallel_for(0UL, num_buckets, [&](size_t i) {
      std::mt19937 local_rng(seeds[i]);
      std::shuffle(shuffled_elements.begin() + bucket_bounds[i], shuffled_elements.begin() + bucket_bounds[i] + 1, local_rng);
    });
  }

  vec<T> elements;
  vec<T> shuffled_elements;

  // convenience
  typename vec<T>::const_iterator begin() const { return shuffled_elements.cbegin(); }
  typename vec<T>::const_iterator end() const { return shuffled_elements.cend(); }
};

template<class T>
class DeterministicParallelUniformRandomPermutation : public BucketSeeding  {
public:
  struct IntegerRange {
    size_t a, b;
    using value_type = size_t;
    size_t operator[](size_t i) const { return a + i;  }
    size_t size() const { return b - a; }
  };

  template<typename F>
  void create_permutation(size_t n, F get_bucket, size_t num_tasks, std::mt19937& rng) {
    permutation.resize(n);

    if (n < 1 << 17) {
      std::iota(permutation.begin(), permutation.end(), 0);
      std::shuffle(permutation.begin(), permutation.end(), rng);
      return;
    }

    vec<uint32_t> bucket_bounds = parallel::counting_sort(IntegerRange{0, n}, permutation, num_buckets, get_bucket, num_tasks);
    ASSERT(bucket_bounds.size() == num_buckets + 1);

    fill_seeds(rng);

    tbb::parallel_for(0UL, num_buckets, [&](size_t i) {
      std::mt19937 local_rng(seeds[i]);
      std::shuffle(permutation.begin() + bucket_bounds[i], permutation.begin() + bucket_bounds[i] + 1, local_rng);
    });
  }

  vec<T> permutation;

  // convenience
  typename vec<T>::const_iterator begin() const { return permutation.cbegin(); }
  typename vec<T>::const_iterator end() const { return permutation.cend(); }
};

struct BucketPrecomputation {

  void compute_buckets(size_t n) {
    // generate more seeds with one starting seed and then generate a bunch of small random numbers from each thread
    precomputed_buckets.resize(n);
  }

  size_t operator()(size_t i) const {
    return precomputed_buckets[i];
  }

  vec<uint8_t> precomputed_buckets;
};

struct BucketHashing {

  void compute_buckets(size_t n, uint32_t seed) {
    unused(n);
    hash.init(seed);
    state = hashing::integer::hash32(seed);
  }

  size_t operator()(uint32_t i) const {
    return hashing::integer::combine(state, hashing::integer::hash32(i));
  }

  hashing::HashTabulated<uint32_t> hash;
  uint32_t state;
};


}
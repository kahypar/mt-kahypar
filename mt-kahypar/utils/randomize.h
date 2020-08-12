/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <limits>
#include <random>
#include <thread>
#include <vector>

#include "tbb/task_group.h"
#include "tbb/parallel_for.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

#include "mt-kahypar/parallel/parallel_counting_sort.h"
#include "hash.h"

namespace mt_kahypar {
namespace utils {

/*!
 * Combines a global seed and an iteration of a loop to initialize an RNG for that iteration
 */
size_t seed_iteration(size_t seed, size_t iteration) {
  return hashing::integer::combine(seed, hashing::integer::hash(iteration));
}

class UniformRandomSelector {
public:
  UniformRandomSelector(size_t hash_function_seed) : rng(hash_function_seed) { }

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
   * Does not reseed the hash function but reseeds the RNG, so that it can be used
   */
  void reset(size_t seed) {
    replace();
    rng.init(seed);
  }

  /*!
   * Use this function to get reproducible pseudorandom numbers in each iteration of a parallel-for loop
   */
  void reset(size_t seed, size_t iteration) {
    reset(seed_iteration(seed, iteration));
  }

private:
  // We want reproducible random numbers with work-stealing. There are two options:
  // Since it is too slow to initialize a std::mt19937 for every vertex, we use a hashing-based RNG
  hashing::SimpleHashRNG<size_t> rng;
  std::uniform_int_distribution<size_t> dist;
  size_t counter = 0;
};

class ParallelSeeding {
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
class DeterministicParallelUniformRandomShuffler : public ParallelSeeding {
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
  vec<T>::const_iterator begin() const { return shuffled_elements.cbegin(); }
  vec<T>::const_iterator end() const { return shuffled_elements.cend(); }
};

template<class T>
class DeterministicParallelUniformRandomPermutation : public ParallelSeeding  {
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
  vec<T>::const_iterator begin() const { return permutation.cbegin(); }
  vec<T>::const_iterator end() const { return permutation.cend(); }
};

class BucketPrecomputation {
public:

  void compute_buckets(size_t n) {
    // generate more seeds with one starting seed and then generate a bunch of small random numbers from each thread
    precomputed_buckets.resize(n);
  }

  size_t operator()(size_t i) const {
    return precomputed_buckets[i];
  }


  vec<uint8_t> precomputed_buckets;
};

class BucketHashing {

  void compute_buckets(size_t n, uint32_t seed) {
    hash.init(seed);
    state = hashing::integer::hash32(seed);
  }

  size_t operator()(uint32_t i) const {
    return hashing::integer::combine(state, hashing::integer::hash32(i));
  }

  hashing::HashTabulated<uint32_t, uint32_t> hash;
  uint32_t state;
};




class Randomize {
  static constexpr bool debug = false;
  static constexpr size_t PRECOMPUTED_FLIP_COINS = 128;

  using SwapBlock = std::pair<size_t, size_t>;

  class RandomFunctions {
   public:
    RandomFunctions() :
      _seed(-1),
      _gen(),
      _next_coin_flip(0),
      _precomputed_flip_coins(PRECOMPUTED_FLIP_COINS),
      _int_dist(0, std::numeric_limits<int>::max()),
      _float_dist(0, 1),
      _norm_dist(0, 1) {
      precompute_flip_coins();
    }

    void setSeed(int seed) {
      _seed = seed;
      _gen.seed(_seed);
      precompute_flip_coins();
    }

    bool flipCoin() {
      return _precomputed_flip_coins[++_next_coin_flip % PRECOMPUTED_FLIP_COINS];
    }

    // returns uniformly random int from the interval [low, high]
    int getRandomInt(int low, int high) {
      return _int_dist(_gen, std::uniform_int_distribution<int>::param_type(low, high));
    }

    // returns uniformly random float from the interval [low, high)
    float getRandomFloat(float low, float high) {
      return _float_dist(_gen, std::uniform_real_distribution<float>::param_type(low, high));
    }

    float getNormalDistributedFloat(float mean, float std_dev) {
      return _norm_dist(_gen, std::normal_distribution<float>::param_type(mean, std_dev));
    }

    std::mt19937& getGenerator() {
      return _gen;
    }

   private:
    void precompute_flip_coins() {
      std::uniform_int_distribution<int> bool_dist(0,1);
      for (size_t i = 0; i < PRECOMPUTED_FLIP_COINS; ++i) {
        _precomputed_flip_coins[i] = static_cast<bool>(bool_dist(_gen));
      }
    }

    int _seed;
    std::mt19937 _gen;
    size_t _next_coin_flip;
    std::vector<bool> _precomputed_flip_coins;
    std::uniform_int_distribution<int> _int_dist;
    std::uniform_real_distribution<float> _float_dist;
    std::normal_distribution<float> _norm_dist;
  };

 public:
  static Randomize& instance() {
    static Randomize instance;
    return instance;
  }

  void enableLocalizedParallelShuffle(const size_t localized_random_shuffle_block_size) {
    _perform_localized_random_shuffle = true;
    _localized_random_shuffle_block_size = localized_random_shuffle_block_size;
  }

  void setSeed(int seed) {
    for (uint32_t i = 0; i < std::thread::hardware_concurrency(); ++i) {
      _rand[i].setSeed(seed + i);
    }
  }

  bool flipCoin(int cpu_id) {
    ASSERT(cpu_id < (int)std::thread::hardware_concurrency());
    return _rand[cpu_id].flipCoin();
  }

  template <typename T>
  void shuffleVector(std::vector<T>& vector, size_t num_elements, int cpu_id) {
    ASSERT(cpu_id < (int)std::thread::hardware_concurrency());
    std::shuffle(vector.begin(), vector.begin() + num_elements, _rand[cpu_id].getGenerator());
  }

  template <typename T>
  void shuffleVector(std::vector<T>& vector, int cpu_id = -1) {
    if (cpu_id == -1)
      cpu_id = sched_getcpu();
    ASSERT(cpu_id < (int)std::thread::hardware_concurrency());
    std::shuffle(vector.begin(), vector.end(), _rand[cpu_id].getGenerator());
  }

  template <typename T>
  void shuffleVector(parallel::scalable_vector<T>& vector, int cpu_id = -1) {
    if (cpu_id == -1)
      cpu_id = sched_getcpu();
    ASSERT(cpu_id < (int)std::thread::hardware_concurrency());
    std::shuffle(vector.begin(), vector.end(), _rand[cpu_id].getGenerator());
  }

  template <typename T>
  void shuffleVector(parallel::scalable_vector<T>& vector, size_t num_elements, int cpu_id) {
    ASSERT(cpu_id < (int)std::thread::hardware_concurrency());
    std::shuffle(vector.begin(), vector.begin() + num_elements, _rand[cpu_id].getGenerator());
  }

  template <typename T>
  void shuffleVector(std::vector<T>& vector, size_t i, size_t j, int cpu_id) {
    ASSERT(i <= j && j <= vector.size());
    ASSERT(cpu_id < (int)std::thread::hardware_concurrency());
    std::shuffle(vector.begin() + i, vector.begin() + j, _rand[cpu_id].getGenerator());
  }

  template <typename T>
  void shuffleVector(parallel::scalable_vector<T>& vector, size_t i, size_t j, int cpu_id) {
    ASSERT(i <= j && j <= vector.size());
    ASSERT(cpu_id < (int)std::thread::hardware_concurrency());
    if ( _perform_localized_random_shuffle ) {
      localizedShuffleVector(vector, i, j, cpu_id);
    } else {
      std::shuffle(vector.begin() + i, vector.begin() + j, _rand[cpu_id].getGenerator());
    }
  }

  template <typename T>
  void parallelShuffleVector(parallel::scalable_vector<T>& vector, const size_t i, const size_t j) {
    ASSERT(i <= j && j <= vector.size());
    const size_t P = 2 * std::thread::hardware_concurrency();
    const size_t N = j - i;
    const size_t step = N / P;

    if ( _perform_localized_random_shuffle ) {
      tbb::parallel_for(0UL, P, [&](const size_t k) {
        const size_t start = i + k * step;
        const size_t end = i + (k == P - 1 ? N : (k + 1) * step);
        localizedShuffleVector(vector, start, end, sched_getcpu());
      });
    } else {
      // Compute blocks that should be swapped before
      // random shuffling
      parallel::scalable_vector<SwapBlock> swap_blocks;
      parallel::scalable_vector<bool> matched_blocks(P, false);
      int cpu_id = sched_getcpu();
      for ( size_t a = 0; a < P; ++a ) {
        if ( !matched_blocks[a] ) {
          matched_blocks[a] = true;
          size_t b = getRandomInt(0, P - 1, cpu_id);
          while ( matched_blocks[b] ) {
            b = ( b + 1 ) % P;
          }
          matched_blocks[b] = true;
          swap_blocks.push_back(std::make_pair(a, b));
        }
      }
      ASSERT(swap_blocks.size() == P / 2, V(swap_blocks.size()) << V(P));

      tbb::parallel_for(0UL, P / 2, [&](const size_t k) {
        const size_t block_1 = swap_blocks[k].first;
        const size_t block_2 = swap_blocks[k].second;
        const size_t start_1 = i + block_1 * step;
        const size_t end_1 = i + (block_1 == P - 1 ? N : (block_1 + 1) * step);
        const size_t start_2 = i + block_2 * step;
        const size_t end_2 = i + (block_2 == P - 1 ? N : (block_2 + 1) * step);
        const int cpu_id = sched_getcpu();
        swapBlocks(vector, start_1, end_1, start_2, end_2);
        std::shuffle(vector.begin() + start_1, vector.begin() + end_1, _rand[cpu_id].getGenerator());
        std::shuffle(vector.begin() + start_2, vector.begin() + end_2, _rand[cpu_id].getGenerator());
      });
    }
  }

  // returns uniformly random int from the interval [low, high]
  int getRandomInt(int low, int high, int cpu_id) {
    ASSERT(cpu_id < (int)std::thread::hardware_concurrency());
    return _rand[cpu_id].getRandomInt(low, high);
  }

  // returns uniformly random float from the interval [low, high)
  float getRandomFloat(float low, float high, int cpu_id) {
    ASSERT(cpu_id < (int)std::thread::hardware_concurrency());
    return _rand[cpu_id].getRandomFloat(low, high);
  }

  float getNormalDistributedFloat(float mean, float std_dev, int cpu_id) {
    ASSERT(cpu_id < (int)std::thread::hardware_concurrency());
    return _rand[cpu_id].getNormalDistributedFloat(mean, std_dev);
  }

  std::mt19937& getGenerator() {
    int cpu_id = sched_getcpu();
    return _rand[cpu_id].getGenerator();
  }

 private:
  explicit Randomize() :
    _rand(std::thread::hardware_concurrency()),
    _perform_localized_random_shuffle(false),
    _localized_random_shuffle_block_size(1024) { }

  template <typename T>
  void swapBlocks(parallel::scalable_vector<T>& vector,
                  const size_t start_1,
                  const size_t end_1,
                  const size_t start_2,
                  const size_t end_2) {
    ASSERT(start_1 <= end_1);
    ASSERT(start_2 <= end_2);
    ASSERT(end_1 <= vector.size());
    ASSERT(end_2 <= vector.size());
    size_t N = std::min(end_1 - start_1, end_2 - start_2);
    for ( size_t i = 0; i < N; ++i ) {
      std::swap(vector[start_1 + i], vector[start_2 + i]);
    }
  }

  template <typename T>
  void localizedShuffleVector(parallel::scalable_vector<T>& vector, const size_t i, const size_t j, const int cpu_id) {
    ASSERT(i <= j && j <= vector.size());
    for ( size_t start = i; start < j; start += _localized_random_shuffle_block_size ) {
      const size_t end = std::min(start + _localized_random_shuffle_block_size, j);
      std::shuffle(vector.begin() + start, vector.begin() + end, _rand[cpu_id].getGenerator());
    }
  }

  std::vector<RandomFunctions> _rand;
  bool _perform_localized_random_shuffle;
  size_t _localized_random_shuffle_block_size;
};
}  // namespace utils
}  // namespace mt_kahypar

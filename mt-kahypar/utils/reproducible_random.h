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
#include <tbb/tick_count.h>

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


struct PrecomputeBucket {
  void compute_buckets(size_t n, size_t num_tasks, uint32_t seed) {
    if (n > precomputed_buckets.size()) {
      precomputed_buckets.resize(n);
    }

    const size_t chunk_size = parallel::chunking::idiv_ceil(n, num_tasks);

    tbb::parallel_for(0UL, num_tasks, [&](size_t i) {
      std::mt19937 rng(seed);
      rng.discard(i);
      rng.seed(rng());

      auto [begin, end] = parallel::chunking::bounds(i, n, chunk_size);
      for (size_t j = begin; j < end; ++j) {
        precomputed_buckets[j] = static_cast<uint8_t>(rng());
      }
    });
  }

  size_t operator()(size_t i) const {
    return precomputed_buckets[i];
  }

  vec<uint8_t> precomputed_buckets;
};

// optimized version of PrecomputeBucket that uses all 32 bits of the rng call
struct PrecomputeBucketOpt {
  void compute_buckets(size_t n, size_t num_tasks, uint32_t seed) {
    if (n > precomputed_buckets.size()) {
      precomputed_buckets.resize(n);
    }

    size_t chunk_size = parallel::chunking::idiv_ceil(n, num_tasks);
    chunk_size += chunk_size % 4; // round up to multiple of 4 --> only last range has to do the overhang bit

    tbb::parallel_for(0UL, num_tasks, [&](size_t i) {
      std::mt19937 rng(seed);
      rng.discard(i);
      rng.seed(rng());

      auto [begin, end] = parallel::chunking::bounds(i, n, chunk_size);

      size_t overhang = end % 4;
      size_t truncated_end = end - overhang;

      for (size_t j = begin; j < truncated_end; j += 4) {
        uint32_t x = rng();
        *( reinterpret_cast<uint32_t*>(precomputed_buckets.data() + j) ) = x;
      }

      if (overhang > 0) {
        uint32_t x = rng();
        for (size_t j = 0; j < overhang; ++j) {
          assert(end - j - 1 < precomputed_buckets.size());
          precomputed_buckets[end - j - 1] = static_cast<uint8_t>( mask & (x >> (8*j)) );
        }
      }
    });
  }

  size_t operator()(size_t i) const {
    return precomputed_buckets[i];
  }

  static constexpr uint32_t mask = (1 << 8) - 1;
  vec<uint8_t> precomputed_buckets;
};


struct BucketHashing {
  void compute_buckets(size_t n, size_t /*num_tasks*/, uint32_t seed) {
    unused(n);
    state = hashing::integer::hash32(seed);
  }

  size_t operator()(uint32_t i) const {
   return hashing::integer::combine32(state, hashing::integer::hash32(i)) % num_buckets;
  }

  static constexpr size_t num_buckets = 256;
  uint32_t state;
};


template<typename T, typename GetBucketCallable = PrecomputeBucketOpt>
class ParallelShuffle {
public:
  vec<T> permutation;
  // convenience
  typename vec<T>::const_iterator begin() const { return permutation.cbegin(); }
  typename vec<T>::const_iterator end() const { return permutation.cend(); }


  template<typename RangeT>
  void shuffle(const RangeT& input_elements, size_t num_tasks, std::mt19937& rng) {
    static_assert(std::is_same<typename RangeT::value_type, T>::value);
    const size_t n = input_elements.size();
    auto t_alloc = tbb::tick_count::now();
    permutation.resize(n);
    if (n < 1 << 15) {
      for (size_t i = 0; i < n; ++i) { permutation[i] = input_elements[i]; }
      std::shuffle(permutation.begin(), permutation.end(), rng);
    } else {
      auto t_comp_buckets = tbb::tick_count::now();
      // assign random buckets to elements. either hash based on position or random tags with static load balancing
      get_bucket.compute_buckets(n, num_tasks, rng());

      auto t_sort = tbb::tick_count::now();
      // sort elements by random buckets
      vec<uint32_t> bucket_bounds = parallel::counting_sort(input_elements, permutation, num_buckets, get_bucket, num_tasks);
      assert(bucket_bounds.size() == num_buckets + 1);

      auto t_shuffle = tbb::tick_count::now();
      // shuffle each bucket
      for (size_t i = 0; i < num_buckets; ++i) {
        seeds[i] = rng();
      }
      tbb::parallel_for(0UL, num_buckets, [&](size_t i) {
        std::mt19937 local_rng(seeds[i]);    // alternative: seed with hash of seed and range begin
        std::shuffle(permutation.begin() + bucket_bounds[i], permutation.begin() + bucket_bounds[i + 1], local_rng);
      });

      auto t_end = tbb::tick_count::now();

      std::cout << "shuffle times: " << (t_end - t_alloc).seconds() << "\n"
                << "\t alloc " << (t_comp_buckets - t_alloc).seconds() << "\n"
                << "\t comp buckets " << (t_sort - t_comp_buckets).seconds() << "\n"
                << "\t sort " << (t_shuffle - t_sort).seconds() << "\n"
                << "\t shuffle " << (t_end - t_shuffle).seconds() << "\n"
                << std::flush;
    }
  }

protected:
  GetBucketCallable get_bucket;
  static constexpr size_t num_buckets = 256;
  std::array<std::mt19937::result_type, num_buckets> seeds;
};


// TODO make prng exchangeable via template parameters everywhere

template<typename IntegralT, typename GetBucketCallable = PrecomputeBucketOpt>
class ParallelPermutation : public ParallelShuffle<IntegralT, GetBucketCallable> {
public:
  void create_integer_permutation(IntegralT n, size_t num_tasks, std::mt19937& rng) {
    static_assert(std::is_integral<IntegralT>::value);
    IntegerRange iota = {0, n};
    this->shuffle(iota, num_tasks, rng);
  }

protected:
  struct IntegerRange {
    IntegralT a, b;
    using value_type = IntegralT;
    IntegralT operator[](IntegralT i) const { return a + i;  }
    size_t size() const { return b - a; }
  };
};

class FeistelPermutation {
public:
  FeistelPermutation(size_t num_rounds, size_t num_entries, std::mt19937& rng) {
    create_permutation(num_rounds, num_entries, rng);
  }

  void create_permutation(size_t num_rounds, size_t num_entries, std::mt19937& rng) {
    // gen keys
    keys.clear();
    for (size_t i = 0; i < num_rounds; ++i) {
      keys.push_back(rng());
    }

    // set bit masks
    uint64_t num_bits = 0, next_power = 1;
    while (next_power < num_entries) {
      next_power <<= 1;
      num_bits++;
    }

    num_right_bits = (num_bits / 2) + (num_bits % 2);
    num_left_bits = num_bits / 2;

    right_mask = (1 << num_right_bits) - 1;
    left_mask = ((1 << num_left_bits) - 1) << num_right_bits;
  }



  uint64_t encrypt(uint64_t x) const {
    uint32_t l = x >> 32;
    uint32_t r = x << 32 >> 32;
    uint32_t next_l;

    for (uint32_t key : keys) {
      next_l = r;
      r = round_function(key, r) ^ l;
      l = next_l;
    }

    // applying mask only to last r and l is bad because we cannot recover those for decryption

    return uint64_t(r) << 32 | uint64_t(l);
  }

  uint64_t decrypt(uint64_t x) const {
    uint32_t l = x >> 32;
    uint32_t r = x << 32 >> 32;
    uint32_t next_l;

    for (auto it = keys.rbegin(); it != keys.rend(); ++it) {
      uint32_t key = *it;
      next_l = r;
      r = round_function(key, r) ^ l;
      l = next_l;
    }

    return uint64_t(r) << 32 | uint64_t(l);
  }

private:

  uint32_t round_function(uint32_t key, uint32_t raw) const {
    return hashing::integer::combine32(key, hashing::integer::hash32(raw));
  }

  vec<uint32_t> keys;
  uint64_t num_right_bits, num_left_bits, right_mask, left_mask;

};

} // namespace mt_kahypar::utils
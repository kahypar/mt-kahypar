/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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
#include <mutex>
#include <random>
#include <thread>
#include <vector>

#include "tbb/task_group.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace utils {
class Randomize {
  static constexpr bool debug = false;
  static constexpr size_t PRECOMPUTED_FLIP_COINS = 128;

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
      bool coinFlip = _precomputed_flip_coins[_next_coin_flip++];
      _next_coin_flip = (_next_coin_flip % PRECOMPUTED_FLIP_COINS);
      return coinFlip;
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
      std::uniform_int_distribution<int> bool_dist;
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
    ASSERT(cpu_id < (int)std::thread::hardware_concurrency());
    std::shuffle(vector.begin() + i, vector.begin() + j, _rand[cpu_id].getGenerator());
  }

  template <typename T>
  void shuffleVector(parallel::scalable_vector<T>& vector, size_t i, size_t j, int cpu_id) {
    ASSERT(cpu_id < (int)std::thread::hardware_concurrency());
    std::shuffle(vector.begin() + i, vector.begin() + j, _rand[cpu_id].getGenerator());
  }

  template <typename T>
  void parallelShuffleVector(parallel::scalable_vector<T>& vector) {
    const size_t P = std::thread::hardware_concurrency();
    const size_t N = vector.size();
    const size_t step = N / P;
    tbb::parallel_for(0UL, P, [&](const size_t i) {
      const size_t start = i * step;
      const size_t end = i == P - 1 ? N : (i + 1) * step;
      const int cpu_id = sched_getcpu();
      std::shuffle(vector.begin() + start, vector.begin() + end, _rand[cpu_id].getGenerator());
    });
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

 private:
  explicit Randomize() :
    _rand(std::thread::hardware_concurrency()) { }

  std::vector<RandomFunctions> _rand;
};
}  // namespace utils
}  // namespace mt_kahypar

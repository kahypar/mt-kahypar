#include "mt-kahypar/utils/reproducible_random.h"

#include <iostream>
#include <random>
#include <tbb/task_scheduler_init.h>

namespace mt_kahypar::utils {


void benchShuffle(size_t n, int num_threads) {
  tbb::task_scheduler_init tsi(num_threads);

#ifndef NDEBUG
  auto is_permutation = [&](vec<int>& r1, vec<int>& r2) {
    std::sort(r1.begin(), r1.end());
    std::sort(r2.begin(), r2.end());
    return r1 == r2;
  };
#endif

  ParallelPermutation<int, PrecomputeBucket> shuffle_preassign;
  ParallelPermutation<int, PrecomputeBucketOpt> shuffle_preassign_opt;
  ParallelPermutation<int, BucketHashing> shuffle_hash;

  uint32_t seed = 420;
  std::mt19937 rng(seed);

  vec<int> comp(n);
  std::iota(comp.begin(), comp.end(), 0);

  std::cout << "preassign buckets" << std::endl;
  shuffle_preassign.create_integer_permutation(n, num_threads, rng);
  assert(is_permutation(comp, shuffle_preassign.permutation));

  std::cout << "preassign buckets opt" << std::endl;
  rng.seed(seed);
  shuffle_preassign_opt.create_integer_permutation(n, num_threads, rng);
  assert(is_permutation(comp, shuffle_preassign_opt.permutation));

  std::cout << "hash" << std::endl;
  rng.seed(seed);
  shuffle_hash.create_integer_permutation(n, num_threads, rng);
  assert(is_permutation(comp, shuffle_hash.permutation));
}



void testFeistel() {
  std::mt19937 rng(420);
  FeistelPermutation feistel(120, rng);

  auto t = [&](uint64_t plain_text) {
    uint64_t encrypted = feistel.encrypt(plain_text);
    uint64_t decrypted = feistel.decrypt(encrypted);
    assert(decrypted == plain_text);
  };

  t(420);
  t(245252);
  t(11);
  t(16841361);
  for (size_t i = 0; i < 500; ++i) {
    t(i);
  }
}


}


int main(int argc, char* argv[]) {
/*
  if (argc != 3) {
    std::cout << "Usage. num-threads permutation-size" << std::endl;
    std::exit(0);
  }

  int num_threads = std::stoi(argv[1]);
  size_t n = std::stoi(argv[2]);
  mt_kahypar::utils::benchShuffle(n, num_threads);
  */
  mt_kahypar::utils::testFeistel();
  return 0;
}

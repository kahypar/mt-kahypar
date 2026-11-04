/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2025 Nikolai Maas <nikolai.maas@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#pragma once

#include <cstring>

#ifdef KAHYPAR_ENABLE_THP
#include <new>
#include <limits>
#include <sys/mman.h>
#endif

#ifndef NDEBUG
#include <memory>
#include <stdlib.h>
#else
#include <tbb/tbb_allocator.h>
#include <tbb/scalable_allocator.h>
#endif

namespace mt_kahypar {
namespace parallel {

#ifndef NDEBUG
// note: tbb allocators will still be used for TBB data structures such as
// tbb::concurrent_vector, tbb::concurrent_queue or tbb::enumerable_thread_specific

template<typename T>
using scalable_allocator = std::allocator<T>;

inline void* mtk_scalable_malloc(size_t size) {
  return malloc(size);
}

inline void* mtk_scalable_calloc(size_t nobj, size_t size) {
  return calloc(nobj, size);
}

inline void mtk_scalable_free(void* ptr) {
  free(ptr);
}

#ifdef KAHYPAR_ENABLE_THP
inline void *mtk_scalable_malloc_thp(size_t size) {
  static constexpr size_t kTHPThreshold = 16 * (1 << 20);
  static constexpr size_t kTHPPageSize = 2 * (1 << 20);

  if (size < kTHPThreshold) {
    return mtk_scalable_malloc(size);
  }

  void *ptr = nullptr;
  posix_memalign(&ptr, kTHPPageSize, size);
  madvise(ptr, size, MADV_HUGEPAGE);

  return ptr;
}
#endif // KAHYPAR_ENABLE_THP
#else // !NDEBUG
template <typename T>
using scalable_allocator = tbb::tbb_allocator<T>;

inline void* mtk_scalable_malloc(size_t size) {
  return scalable_malloc(size);
}

inline void* mtk_scalable_calloc(size_t nobj, size_t size) {
  return scalable_calloc(nobj, size);
}

inline void mtk_scalable_free(void* ptr) {
  scalable_free(ptr);
}
#ifdef KAHYPAR_ENABLE_THP
inline void* mtk_scalable_malloc_thp(size_t size) {
  static constexpr size_t kTHPThreshold = 16 * (1 << 20);
  static constexpr size_t kTHPPageSize = 2 * (1 << 20);

  if (size < kTHPThreshold) {
    return mtk_scalable_malloc(size);
  }

  void* ptr = nullptr;
  scalable_posix_memalign(&ptr, kTHPPageSize, size);
  madvise(ptr, size, MADV_HUGEPAGE);

  return ptr;
}
#endif // KAHYPAR_ENABLE_THP
#endif // !NDEBUG
#ifdef KAHYPAR_ENABLE_THP
inline void *mtk_scalable_calloc_thp(size_t nobj, size_t size) {
  const size_t nbytes = nobj * size;

  void *ptr = mtk_scalable_malloc_thp(nbytes);
  if (nullptr == ptr) {
    return ptr;
  }
  std::memset(ptr, 0, nbytes);

  return ptr;
}
#else // KAHYPAR_ENABLE_THP
inline void* mtk_scalable_malloc_thp(size_t size) {
  return mtk_scalable_malloc(size);
}
inline void* mtk_scalable_calloc_thp(size_t nobj, size_t size) {
  return mtk_scalable_calloc(nobj, size);
}
#endif // KAHYPAR_ENABLE_THP

inline void mtk_scalable_free_thp(void *ptr) {
  mtk_scalable_free(ptr);
}

#ifdef KAHYPAR_ENABLE_THP
template <typename T, bool zero_initialize = false>
class thp_allocator {
public:
  using value_type = T;

  thp_allocator() noexcept = default;

  template <typename U> thp_allocator(const thp_allocator<U> &) noexcept {}

  [[nodiscard]] T *allocate(std::size_t n) {
    if (n > std::numeric_limits<std::size_t>::max() / sizeof(T)) {
      throw std::bad_array_new_length();
    }

    size_t nbytes = n * sizeof(T);

    if (nbytes < kTHPThreshold) {
      return scalable_allocator<T>{}.allocate(n);
    }

    T *ptr = nullptr;
    if constexpr (zero_initialize) {
      ptr = mtk_scalable_calloc_thp(1, nbytes);
    } else {
      ptr = mtk_scalable_malloc_thp(nbytes);
    }
    if (nullptr == ptr) {
      throw std::bad_alloc();
    }

    return static_cast<T *>(ptr);
  }

  void deallocate(T *p, std::size_t n) noexcept {
    size_t size = n * sizeof(T);

    if (size < kTHPThreshold) {
      scalable_allocator<T>{}.deallocate(p, n);
      return;
    }

    mtk_scalable_free_thp(p);
  }

private:
  static constexpr size_t kTHPThreshold = 16 * (1 << 20);
  static constexpr size_t kTHPPageSize = 2 * (1 << 20);
};

template <typename T>
using scalable_thp_allocator = thp_allocator<T>;
#else // KAHYPAR_ENABLE_THP
template<typename T>
using scalable_thp_allocator = scalable_allocator<T>;
#endif // KAHYPAR_ENABLE_THP

}  // namespace parallel
}  // namespace mt_kahypar

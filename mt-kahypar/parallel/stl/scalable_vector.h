/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <type_traits>
#include <vector>

#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>

#include "mt-kahypar/parallel/stl/allocator.h"
#include "mt-kahypar/macros.h"

namespace mt_kahypar {

template<typename T>
using vec = std::vector<T, parallel::scalable_allocator<T> >;  // shorter name

namespace parallel {
template <typename T>
using scalable_vector = vec<T>;

template <typename>
struct is_scalable_vector : std::false_type {};
template <typename T>
struct is_scalable_vector<scalable_vector<T>> : std::true_type {};
#ifdef KAHYPAR_ENABLE_THP
template <typename T>
struct is_scalable_vector<scalable_thp_vector<T>> : std::true_type {};
#endif

template <typename T>
concept ScalableVector = is_scalable_vector<std::remove_cvref_t<T>>::value;

template<ScalableVector T>
static inline void free(T& vec) {
  T tmp_vec;
  vec = std::move(tmp_vec);
}

template<ScalableVector T> requires ScalableVector<typename T::value_type>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static void parallel_free(T& vec) {
  tbb::parallel_for(UL(0), vec.size(), [&](const size_t i) {
    free(vec[i]);
  });
}

template<ScalableVector... ScalableVectors>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static void parallel_free(ScalableVectors&... vecs) {
    tbb::parallel_invoke(
        [&] { free(vecs); }...
    );
}

}  // namespace parallel


}  // namespace mt_kahypar

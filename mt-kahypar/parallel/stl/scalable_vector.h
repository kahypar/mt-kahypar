/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <vector>

#include "tbb/parallel_for.h"
#include "tbb/parallel_invoke.h"
#include "tbb/scalable_allocator.h"

#include "mt-kahypar/macros.h"

namespace mt_kahypar {

template<typename T>
using vec = std::vector<T, tbb::scalable_allocator<T> >;  // shorter name

namespace parallel {
template <typename T>
using scalable_vector = std::vector<T, tbb::scalable_allocator<T> >;

template<typename T>
static inline void free(scalable_vector<T>& vec) {
  scalable_vector<T> tmp_vec;
  vec = std::move(tmp_vec);
}

template<typename T>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static void parallel_free(scalable_vector<scalable_vector<T>>& vec) {
  tbb::parallel_for(0UL, vec.size(), [&](const size_t i) {
    free(vec[i]);
  });
}

template<typename T1,
         typename T2>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static void parallel_free(scalable_vector<T1>& vec1,
                                                             scalable_vector<T2>& vec2) {
  tbb::parallel_invoke([&] {
    free(vec1);
  }, [&] {
    free(vec2);
  });
}

template<typename T1,
         typename T2,
         typename T3>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static void parallel_free(scalable_vector<T1>& vec1,
                                                             scalable_vector<T2>& vec2,
                                                             scalable_vector<T3>& vec3) {
  tbb::parallel_invoke([&] {
    free(vec1);
  }, [&] {
    free(vec2);
  }, [&] {
    free(vec3);
  });
}

template<typename T1,
         typename T2,
         typename T3,
         typename T4>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static void parallel_free(scalable_vector<T1>& vec1,
                                                             scalable_vector<T2>& vec2,
                                                             scalable_vector<T3>& vec3,
                                                             scalable_vector<T4>& vec4) {
  tbb::parallel_invoke([&] {
    free(vec1);
  }, [&] {
    free(vec2);
  }, [&] {
    free(vec3);
  }, [&] {
    free(vec4);
  });
}


template<typename T1,
         typename T2,
         typename T3,
         typename T4,
         typename T5>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static void parallel_free(scalable_vector<T1>& vec1,
                                                             scalable_vector<T2>& vec2,
                                                             scalable_vector<T3>& vec3,
                                                             scalable_vector<T4>& vec4,
                                                             scalable_vector<T5>& vec5) {
  tbb::parallel_invoke([&] {
    free(vec1);
  }, [&] {
    free(vec2);
  }, [&] {
    free(vec3);
  }, [&] {
    free(vec4);
  }, [&] {
    free(vec5);
  });
}

template<typename T1,
         typename T2,
         typename T3,
         typename T4,
         typename T5,
         typename T6>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static void parallel_free(scalable_vector<T1>& vec1,
                                                             scalable_vector<T2>& vec2,
                                                             scalable_vector<T3>& vec3,
                                                             scalable_vector<T4>& vec4,
                                                             scalable_vector<T5>& vec5,
                                                             scalable_vector<T6>& vec6) {
  tbb::parallel_invoke([&] {
    free(vec1);
  }, [&] {
    free(vec2);
  }, [&] {
    free(vec3);
  }, [&] {
    free(vec4);
  }, [&] {
    free(vec5);
  }, [&] {
    free(vec6);
  });
}

}  // namespace parallel


}  // namespace mt_kahypar

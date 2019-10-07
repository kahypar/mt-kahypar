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

#include "kahypar/datastructure/fast_reset_flag_array.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace ds {

template <typename UnderlyingType = std::uint16_t>
class ThreadSafeFastResetFlagArray {

  static constexpr bool debug = false;

 public:
  ThreadSafeFastResetFlagArray(const size_t num_cpus,
                               const size_t size) :
    _data() {
    for ( size_t i = 0; i < num_cpus; ++i ) {
      _data.emplace_back(size);
    }
  }

  ThreadSafeFastResetFlagArray(const ThreadSafeFastResetFlagArray&) = delete;
  ThreadSafeFastResetFlagArray& operator= (const ThreadSafeFastResetFlagArray&) = delete;

  ThreadSafeFastResetFlagArray(ThreadSafeFastResetFlagArray&& other) = default;
  ThreadSafeFastResetFlagArray& operator= (ThreadSafeFastResetFlagArray&&) = default;

  ~ThreadSafeFastResetFlagArray() = default;

  bool get(const size_t cpu_id, const size_t i) const {
    ASSERT(cpu_id < _data.size());
    return _data[cpu_id][i];
  }

  void set(const size_t cpu_id, const size_t i, const bool value) {
    ASSERT(cpu_id < _data.size());
    _data[cpu_id].set(i, value);
  }

  void reset(const size_t cpu_id) {
    ASSERT(cpu_id < _data.size());
    _data[cpu_id].reset();
  }

 private:
  parallel::scalable_vector<kahypar::ds::FastResetFlagArray<UnderlyingType>> _data;

};

} // namespace ds
} // namespace mt_kahypar
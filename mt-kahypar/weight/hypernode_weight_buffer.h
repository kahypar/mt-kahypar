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

#include <memory>

#include <tbb/enumerable_thread_specific.h>

#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/weight/hypernode_weight_base.h"

namespace mt_kahypar {
namespace weight {

class HypernodeWeightBuffer {
  static constexpr size_t CHUNK_SIZE = 8;

  struct Chunk {
    Chunk(Dimension dimension) {
      _data = std::make_unique<weight::HNWeightScalar[]>(dimension * CHUNK_SIZE);
      ENABLE_IF_DEBUG_WEIGHTS(
        _ref_count = std::make_unique<uint32_t[]>(CHUNK_SIZE);
      );
    }

    std::pair<weight::HNWeightScalar*, uint32_t*> allocateAtOffset(Dimension dimension, size_t offset) {
      ASSERT(offset < CHUNK_SIZE);
      weight::HNWeightScalar* data = _data.get() + (dimension * offset);
      uint32_t* ref_count = nullptr;
      ENABLE_IF_DEBUG_WEIGHTS(
        ref_count = _ref_count.get() + offset;
      );
      return {data, ref_count};
    }

    std::unique_ptr<weight::HNWeightScalar[]> _data;

   #ifdef MT_KAHYPAR_DEBUG_WEIGHTS
    std::unique_ptr<uint32_t[]> _ref_count;
   #endif
  };

 public:
  HypernodeWeightBuffer() :
    _local_buffer(),
    _dimension(0) { }

  explicit HypernodeWeightBuffer(const Dimension dimension) :
    _local_buffer(),
    _dimension(dimension) { }

  HypernodeWeightBuffer(const HypernodeWeightBuffer&) = delete;
  HypernodeWeightBuffer & operator= (const HypernodeWeightBuffer&) = delete;

  HypernodeWeightBuffer(HypernodeWeightBuffer&& other) = default;
  HypernodeWeightBuffer& operator=(HypernodeWeightBuffer&& other) = default;

  std::pair<weight::HNWeightScalar*, uint32_t*> allocateAtPosition(size_t pos) {
    ASSERT(_dimension > 0);
    size_t chunk_index = pos / CHUNK_SIZE;
    size_t offset = pos % CHUNK_SIZE;
    auto& buffer = _local_buffer.local();
    while (chunk_index >= buffer.size()) {
      buffer.emplace_back(_dimension);
    }
    return buffer[chunk_index].allocateAtOffset(_dimension, offset);
  }

  Dimension dimension() const {
    return _dimension;
  }

  HypernodeWeightBuffer& getHypernodeWeightBuffer() {
    ASSERT(_dimension > 0);
    return *this;
  }

 private:
  tbb::enumerable_thread_specific<vec<Chunk>> _local_buffer;
  Dimension _dimension;
};


// ###############  BUFFER ALLOCATION FUNCTIONS  ###############

template <typename Proxy>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightRef allocateAt(Proxy& proxy, size_t pos) {
  auto& buffer = proxy.getHypernodeWeightBuffer();
  const Dimension dimension = buffer.dimension();
  const auto [data, ref_count] = buffer.allocateAtPosition(pos);
  ENABLE_IF_DEBUG_WEIGHTS(
    ALWAYS_ASSERT(ref_count != nullptr);
    ALWAYS_ASSERT(*ref_count == 0, "Double allocation at position" << pos);
    // poison the data with garbage values since the content is undefined
    for (Dimension i = 0; i < dimension; i++) {
      data[i] = (std::numeric_limits<HNWeightScalar>::max() / 2) + i;
    }
  )
  return HNWeightRef(data, dimension, ref_count);
}

template <typename Proxy, typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightRef allocateCopy(Proxy& proxy, size_t pos, const Expr& expr) {
  HNWeightRef result = allocateAt(proxy, pos);
  result = expr;
  return result;
}

}  // namespace weight
}  // namespace mt_kahypar

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

#include "mt-kahypar/macros.h"
#include "mt-kahypar/weight/hypernode_weight_base.h"
#include "mt-kahypar/weight/hypernode_weight_operators.h"

namespace mt_kahypar {
namespace weight {

// a single hypernode weight using a separate allocation
class AllocatedHNWeight {
 public:
  AllocatedHNWeight() :
    _data(nullptr),
    _dimension(0) { }

  explicit AllocatedHNWeight(const Dimension dimension,
                             const weight::HNWeightScalar init_value) :
    _data(nullptr),
    _dimension(0) {
      ASSERT(dimension > 0);
      setDimension(dimension, init_value);
    }

  AllocatedHNWeight(const AllocatedHNWeight& other) = delete;

  // assignment operator does complete reassignment
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE AllocatedHNWeight& operator=(const AllocatedHNWeight& other) {
    assign(other);
    return *this;
  }

  AllocatedHNWeight(AllocatedHNWeight&& other) = default;
  AllocatedHNWeight& operator=(AllocatedHNWeight&& other) = default;

  template <typename Expr>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE AllocatedHNWeight& operator=(const Expr& expression) {
    assign(expression);
    return *this;
  }

  template <typename Expr>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE AllocatedHNWeight& operator+=(const Expr& expression) {
    static_cast<HNWeightRef>(*this) += expression;
    return *this;
  }

  template <typename Expr>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE AllocatedHNWeight& operator-=(const Expr& expression) {
    static_cast<HNWeightRef>(*this) -= expression;
    return *this;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE AllocatedHNWeight& operator*=(HNWeightScalar factor) {
    HNWeightRef ref = *this;
    ref *= factor;
    return *this;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE AllocatedHNWeight& operator*=(double factor) {
    HNWeightRef ref = *this;
    ref *= factor;
    return *this;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE AllocatedHNWeight& operator/=(HNWeightScalar divisor) {
    HNWeightRef ref = *this;
    ref /= divisor;
    return *this;
  }

  template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE auto fetch_add(const Expr& right, std::memory_order memory_order) {
    return static_cast<HNWeightRef>(*this).fetch_add(right, memory_order);
  }

  template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE auto add_fetch(const Expr& right, std::memory_order memory_order) {
    return static_cast<HNWeightRef>(*this).add_fetch(right, memory_order);
  }

  template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE auto fetch_sub(const Expr& right, std::memory_order memory_order) {
    return static_cast<HNWeightRef>(*this).fetch_sub(right, memory_order);
  }

  template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE auto sub_fetch(const Expr& right, std::memory_order memory_order) {
    return static_cast<HNWeightRef>(*this).sub_fetch(right, memory_order);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightAtomicCRef load(std::memory_order memory_order) const {
    return static_cast<HNWeightConstRef>(*this).load(memory_order);
  }

  // ####################### Basic Operators #######################

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar at(Dimension i) const {
    ASSERT(i < _dimension);
    return _data[i];
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar& at(Dimension i) {
    ASSERT(i < _dimension);
    return _data[i];
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Dimension dimension() const {
    return _dimension;
  }

  // ####################### Initialization #######################

  void setDimension(const Dimension new_dimension,
                    const weight::HNWeightScalar init_value = 0) {
    ASSERT(_data == nullptr && _dimension == 0 && new_dimension > 0);
    _data = std::make_unique<weight::HNWeightScalar[]>(new_dimension);
    for (size_t i = 0; i < new_dimension; ++i) {
      _data[i] = init_value;
    }
    _dimension = new_dimension;
  }

  // ####################### Conversion #######################

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE operator HNWeightRef() {
    return HNWeightRef(_data.get(), _dimension);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE operator HNWeightConstRef() const {
    return HNWeightConstRef(_data.get(), _dimension);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightRef get() {
    return static_cast<HNWeightRef>(*this);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightConstRef get() const {
    return static_cast<HNWeightConstRef>(*this);
  }

  // ####################### COPY #######################

  AllocatedHNWeight copy() const {
    if (_dimension == 0) return AllocatedHNWeight();

    AllocatedHNWeight result(_dimension, 0);
    result = static_cast<HNWeightConstRef>(*this);
    return result;
  }

 protected:
  template <typename Expr>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void assign(const Expr& expression) {
    if (_data == nullptr) {
      setDimension(expression.dimension());
    }
    static_cast<HNWeightRef>(*this) = expression;
  }

 private:
  std::unique_ptr<weight::HNWeightScalar[]> _data;
  Dimension _dimension;
};
MT_KAHYPAR_DEFINE_VALID_EXPRESSION(AllocatedHNWeight);


class CopyableAllocatedHNWeight: public AllocatedHNWeight {
 public:
  CopyableAllocatedHNWeight() : AllocatedHNWeight() { }

  CopyableAllocatedHNWeight(const Dimension dimension,
                            const weight::HNWeightScalar init_value = 0) :
    AllocatedHNWeight(dimension, init_value) { }

  CopyableAllocatedHNWeight(const CopyableAllocatedHNWeight& other) : AllocatedHNWeight(other.copy()) { }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE CopyableAllocatedHNWeight& operator=(const CopyableAllocatedHNWeight& other) {
    if (dimension() == other.dimension()) {
      assign(other);
    } else {
      static_cast<AllocatedHNWeight&>(*this) = other.copy();
    }
    return *this;
  }

  template <typename Expr>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE CopyableAllocatedHNWeight& operator=(const Expr& expression) {
    assign(expression);
    return *this;
  }

  CopyableAllocatedHNWeight(CopyableAllocatedHNWeight&& other) = default;
  CopyableAllocatedHNWeight& operator=(CopyableAllocatedHNWeight&& other) = default;
};
MT_KAHYPAR_DEFINE_VALID_EXPRESSION(CopyableAllocatedHNWeight);

}  // namespace weight
}  // namespace mt_kahypar

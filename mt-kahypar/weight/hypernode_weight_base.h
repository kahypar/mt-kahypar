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

#include <atomic>
#include <cstdint>
#include <type_traits>
#include <memory>

#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"

namespace mt_kahypar {
namespace weight {

// basic macro definitions
#ifdef KAHYPAR_USE_ASSERTIONS
#define MT_KAHYPAR_DEBUG_WEIGHTS true
#define ENABLE_IF_DEBUG_WEIGHTS(X) X
#else
#define MT_KAHYPAR_DEBUG_WEIGHTS false
#define ENABLE_IF_DEBUG_WEIGHTS(X)
#endif

// basic type definitions
#if KAHYPAR_USE_64_BIT_WEIGHTS
using HNWeightScalar = int64_t;
#else
using HNWeightScalar = int32_t;
#endif
using Dimension = uint32_t;

template<typename T>
struct IsHypernodeweightExpression {
  static constexpr bool value = false;
};

template<typename T>
static constexpr bool is_hypernodeweight_expression = IsHypernodeweightExpression<T>::value;


// helper functions and macros
#define REQUIRE_VALID_WEIGHT(X) typename = std::enable_if_t<weight::is_hypernodeweight_expression<X>>

#define REQUIRE_VALID_WEIGHT_2(X, Y) typename = std::enable_if_t<weight::is_hypernodeweight_expression<X> && weight::is_hypernodeweight_expression<Y>>

#define DEDUCE_TYPE(X) std::remove_reference_t<decltype(X)>

#define MT_KAHYPAR_DEFINE_VALID_EXPRESSION(X) \
  template<> \
  struct IsHypernodeweightExpression<X> { \
    static constexpr bool value = true; \
  }

#define MT_KAHYPAR_DEFINE_VALID_EXPRESSION_UNARY(X) \
  template<typename A> \
  struct IsHypernodeweightExpression<X<A>> { \
    static constexpr bool value = is_hypernodeweight_expression<A>; \
  }

#define MT_KAHYPAR_DEFINE_VALID_EXPRESSION_BINARY(X) \
  template<typename A, typename B> \
  struct IsHypernodeweightExpression<X<A,B>> { \
    static constexpr bool value = is_hypernodeweight_expression<A> && is_hypernodeweight_expression<B>; \
  }

#define CREATE_ATOMIC_HNWEIGHTREF_EXPRESSION(A, B) \
template <typename R> \
class [[nodiscard]] A { \
 public: \
  static_assert(is_hypernodeweight_expression<R>); \
 \
  explicit A(const HNWeightRef& left, const R& right, std::memory_order memory_order): \
    A(left._data, right, memory_order) { \
    ASSERT(left.dimension() == right.dimension()); \
  } \
  \
  A(A&& other): \
    _data(other._data), _right(std::move(other._right)), _memory_order(other._memory_order) { \
    ENABLE_IF_DEBUG_WEIGHTS( \
      _application_tracker = std::move(other._application_tracker); \
      other._application_tracker = nullptr; \
    ) \
  } \
 \
  ENABLE_IF_DEBUG_WEIGHTS( \
    ~A() { \
      if (_application_tracker) { \
        for (Dimension i = 0; i < dimension(); ++i) { \
          ASSERT(_application_tracker[i] == 1, \
            "Atomic expression was not applied once, but:" << _application_tracker[i] << "times; i =" << i); \
        } \
      } \
    } \
  ) \
 \
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar at(Dimension i) const { \
    ENABLE_IF_DEBUG_WEIGHTS( \
      _application_tracker[i] += 1; \
    ) \
    return std::atomic_ref(_data[i]).B(_right.at(i), _memory_order); \
  } \
 \
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Dimension dimension() const { \
    return _right.dimension(); \
  } \
 \
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE A<R> get() const { \
    if constexpr (MT_KAHYPAR_DEBUG_WEIGHTS) { \
      auto result = std::move(const_cast<A<R>&>(*this)); \
      return result; \
    } else { \
      return A<R>(_data, _right, _memory_order); \
    } \
  } \
 \
 private: \
  explicit A(HNWeightScalar* data, const R& right, std::memory_order memory_order): \
    _data(data), _right(right.get()), _memory_order(memory_order) { \
    ENABLE_IF_DEBUG_WEIGHTS( \
      _application_tracker = std::make_shared<int[]>(right.dimension(), 0); \
    ) \
  } \
  \
  HNWeightScalar* _data; \
  R _right; \
  std::memory_order _memory_order; \
  ENABLE_IF_DEBUG_WEIGHTS(std::shared_ptr<int[]> _application_tracker;) \
}


// #################  REFERENCE BASE  #################

template<typename Subclass>
class HNWeightRefBase {
 public:
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Dimension dimension() const {
    return _dimension;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar* get_raw_data() const {
    return _data;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void set_raw_data(HNWeightScalar* data, Dimension dimension) {
    _data = data;
    _dimension = dimension;
  }

 protected:
  template<typename Other>
  friend class HNWeightRefBase;

  explicit HNWeightRefBase(): _data(nullptr), _dimension(0) {}

  explicit HNWeightRefBase(HNWeightScalar* data, Dimension dimension):
    _data(data), _dimension(dimension) {}

  HNWeightRefBase(const HNWeightRefBase& other):
    _data(other._data),
    _dimension(other._dimension) {}

  HNWeightRefBase& operator=(const HNWeightRefBase& other) {
    assignOther(other);
    return *this;
  }

  HNWeightRefBase(HNWeightRefBase&& other):
    _data(other._data),
    _dimension(other._dimension) {
    ENABLE_IF_DEBUG_WEIGHTS(
      other.setEmpty();
    )
  }

  HNWeightRefBase& operator=(HNWeightRefBase&& other) {
    assignOther(other);
    return *this;
  }

  template<typename Other>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void assignOther(const HNWeightRefBase<Other>& other) {
    _data = other._data;
    _dimension = other._dimension;
  }

  template<typename Other>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void assignOther(HNWeightRefBase<Other>&& other) {
    _data = other._data;
    _dimension = other._dimension;
    ENABLE_IF_DEBUG_WEIGHTS(
      other.setEmpty();
    )
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar at(Dimension i) const {
    ASSERT(i < _dimension);
    return _data[i];
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void set(Dimension i, HNWeightScalar value) {
    ASSERT(i < _dimension);
    _data[i] = value;
  }

ENABLE_IF_DEBUG_WEIGHTS(
  void setEmpty() {
    _data = nullptr;
    _dimension = 0;
  }
)

 protected:
  HNWeightScalar* _data;
  Dimension _dimension;
};


// ####################  ATOMIC REFERENCE  ####################

class HNWeightAtomicCRef: public HNWeightRefBase<HNWeightAtomicCRef> {
  using Base = HNWeightRefBase<HNWeightAtomicCRef>;
  using Base::_data;

 public:
  // construction from references
  template<typename Other>
  explicit HNWeightAtomicCRef(HNWeightRefBase<Other>&& other, std::memory_order memory_order):
    Base(), _memory_order(memory_order) {
      assignOther(std::move(other));
    }

  // we need copy construction/assignment for atomics so that the get() contract is fulfilled
  HNWeightAtomicCRef(const HNWeightAtomicCRef&) = default;
  HNWeightAtomicCRef& operator=(const HNWeightAtomicCRef&) = default;

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar at(Dimension i) const {
    ASSERT(i < _dimension);
    return parallel::atomic_load(_data[i], _memory_order);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightAtomicCRef get() const {
    return HNWeightAtomicCRef(*this, _memory_order);
  }

 private:
  template<typename Other>
  explicit HNWeightAtomicCRef(const HNWeightRefBase<Other>& other, std::memory_order memory_order):
    Base(other), _memory_order(memory_order) {}

  std::memory_order _memory_order;
};
MT_KAHYPAR_DEFINE_VALID_EXPRESSION(HNWeightAtomicCRef);


// ####################  CONST REFERENCE  ####################

// forward declaration
class HNWeightRef;

class HNWeightConstRef: public HNWeightRefBase<HNWeightConstRef> {
  using Base = HNWeightRefBase<HNWeightConstRef>;

 public:
  explicit HNWeightConstRef(): Base() {}

  explicit HNWeightConstRef(const HNWeightScalar* data, Dimension dimension):
    Base(const_cast<HNWeightScalar*>(data), dimension) {}

  // copy construction/assignment
  HNWeightConstRef(const HNWeightConstRef&) = default;
  HNWeightConstRef& operator=(const HNWeightConstRef&) = default;

  // construction/assignment from non-const reference
  HNWeightConstRef(const HNWeightRefBase<HNWeightRef>& other): Base() {
    assignOther(std::move(other));
  }

  HNWeightConstRef& operator=(const HNWeightRefBase<HNWeightRef>& other) {
    assignOther(other);
    return *this;
  };

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar at(Dimension i) const {
    return Base::at(i);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightConstRef get() const {
    return HNWeightConstRef(_data, _dimension);  // intentionally no ref count
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightAtomicCRef load(std::memory_order memory_order) const {
    HNWeightConstRef copy = *this;
    return HNWeightAtomicCRef(std::move(copy), memory_order);
  }
};
MT_KAHYPAR_DEFINE_VALID_EXPRESSION(HNWeightConstRef);


// ####################  NON-CONST REFERENCE  ####################

class HNWeightRef: public HNWeightRefBase<HNWeightRef> {
  using Base = HNWeightRefBase<HNWeightRef>;

  // define structs for atomic expressions
  CREATE_ATOMIC_HNWEIGHTREF_EXPRESSION(FetchAddExpr, fetch_add);
  CREATE_ATOMIC_HNWEIGHTREF_EXPRESSION(FetchSubExpr, fetch_sub);

 public:
  explicit HNWeightRef() : Base() {}

  explicit HNWeightRef(HNWeightScalar* data, Dimension dimension):
    Base(data, dimension) {}

  // copy construction is forbidden to avoid aliasing problems
  HNWeightRef(const HNWeightRef&) = delete;

  // move construction is allowed
  HNWeightRef(HNWeightRef&&) = default;

  // assignment operators assign the content, not the reference!
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightRef& operator=(const HNWeightRef& other) {
    return assign(other);
  }

  // assignment operators assign the content, not the reference!
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightRef& operator=(HNWeightRef&& other) {
    return assign(other);
  }

  template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightRef& operator=(const Expr& expression) {
    return assign(expression);
  }

  template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightRef& operator+=(const Expr& expression) {
    return *this = (*this + expression);
  }

  template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightRef& operator-=(const Expr& expression) {
    return *this = (*this - expression);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar at(Dimension i) const {
    return Base::at(i);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar& at(Dimension i) {
    ASSERT(i < _dimension);
    return _data[i];
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void set(Dimension i, HNWeightScalar value) {
    Base::set(i, value);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightConstRef get() const {
    return HNWeightConstRef(_data, _dimension);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightAtomicCRef load(std::memory_order memory_order) && {
    return HNWeightAtomicCRef(std::move(*this), memory_order);
  }

  template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void store(const Expr& right, std::memory_order memory_order) {
    ASSERT(dimension() == right.dimension());
    for (Dimension i = 0; i < dimension(); i++) {
      std::atomic_ref(_data[i]).store(right.at(i), memory_order);
    }
  }

  template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE auto fetch_add(const Expr& right, std::memory_order memory_order) {
    return FetchAddExpr<DEDUCE_TYPE(right.get())>(*this, right.get(), memory_order);
  }

  template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE auto fetch_sub(const Expr& right, std::memory_order memory_order) {
    return FetchSubExpr<DEDUCE_TYPE(right.get())>(*this, right.get(), memory_order);
  }

 private:
  template <typename Expr>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightRef& assign(const Expr& expression) {
    static_assert(is_hypernodeweight_expression<Expr>);
    ASSERT(dimension() == expression.dimension());
    for (Dimension i = 0; i < dimension(); i++) {
      set(i, expression.at(i));
    }
    return *this;
  }
};
MT_KAHYPAR_DEFINE_VALID_EXPRESSION_UNARY(HNWeightRef::FetchAddExpr);
MT_KAHYPAR_DEFINE_VALID_EXPRESSION_UNARY(HNWeightRef::FetchSubExpr);
MT_KAHYPAR_DEFINE_VALID_EXPRESSION(HNWeightRef);

} // namespace weight
} // namespace mt_kahypar

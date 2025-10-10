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

#include <cstdint>
#include <limits>
#include <type_traits>
#include <memory>

#include "mt-kahypar/macros.h"

#ifdef KAHYPAR_USE_ASSERTIONS
#define MT_KAHYPAR_DEBUG_WEIGHTS
#endif

namespace mt_kahypar {

namespace weight {

template<typename T>
struct IsHypernodeweightExpression {
  static constexpr bool value = false;
};

template<typename T>
static constexpr bool is_hypernodeweight_expression = IsHypernodeweightExpression<T>::value;

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE int transformMemoryOrder(std::memory_order memory_order) {
  switch (memory_order) {
    case std::memory_order_relaxed: return __ATOMIC_RELAXED;
    case std::memory_order_consume: return __ATOMIC_CONSUME;
    case std::memory_order_acquire: return __ATOMIC_ACQUIRE;
    case std::memory_order_release: return __ATOMIC_RELEASE;
    case std::memory_order_acq_rel: return __ATOMIC_ACQ_REL;
    case std::memory_order_seq_cst: return __ATOMIC_SEQ_CST;
  };
  return __ATOMIC_SEQ_CST;
}

#ifdef MT_KAHYPAR_DEBUG_WEIGHTS
#define ENABLE_IF_DEBUG_WEIGHTS(X) X
#else
#define ENABLE_IF_DEBUG_WEIGHTS(X)
#endif

#define REQUIRE_VALID_WEIGHT(X) typename = std::enable_if_t<weight::is_hypernodeweight_expression<X>>

#define REQUIRE_VALID_WEIGHT_2(X, Y) typename = std::enable_if_t<weight::is_hypernodeweight_expression<X> && weight::is_hypernodeweight_expression<Y>>

#define DEDUCE_TYPE(X) std::remove_reference_t<decltype(X)>

using HNWeightScalar = int32_t;
using Dimension = uint32_t;

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
  explicit A(const HNWeightRef& left, R&& right, std::memory_order memory_order): \
    _data(left._data), _right(std::move(right)), _memory_order(memory_order) { \
    ASSERT(left.dimension() == right.dimension()); \
    ENABLE_IF_DEBUG_WEIGHTS( \
      _application_tracker = std::make_shared<std::vector<int>>(left.dimension(), 0); \
    ) \
  } \
 \
  ENABLE_IF_DEBUG_WEIGHTS( \
    ~A() { \
      for (Dimension i = 0; i < dimension(); ++i) { \
        ASSERT((*_application_tracker).at(i) == 1, \
          "Atomic expression was not applied once, but:" << (*_application_tracker).at(i)); \
      } \
    } \
  ) \
 \
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar at(Dimension i) const { \
    ENABLE_IF_DEBUG_WEIGHTS( \
      (*_application_tracker).at(i) += 1; \
    ) \
    return B(&_data[i], _right.at(i), transformMemoryOrder(_memory_order)); \
  } \
 \
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Dimension dimension() const { \
    return _right.dimension(); \
  } \
 \
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE A<R> get() const { \
    return *this; \
  } \
 private: \
  HNWeightScalar* _data; \
  R _right; \
  std::memory_order _memory_order; \
  ENABLE_IF_DEBUG_WEIGHTS(std::shared_ptr<std::vector<int>> _application_tracker;) \
}

#define CREATE_BINARY_HNWEIGHT_OPERATOR(A, OPNAME, FUNC) \
template <typename L, typename R> \
class [[nodiscard]] A { \
 public: \
  static_assert(is_hypernodeweight_expression<L> && is_hypernodeweight_expression<R>); \
 \
  explicit A(L&& left, R&& right): _left(std::move(left)), _right(std::move(right)) { \
    ASSERT(left.dimension() == right.dimension()); \
  } \
 \
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar at(Dimension i) const { \
    return FUNC (_left.at(i), _right.at(i)); \
  } \
 \
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Dimension dimension() const { \
    return _left.dimension(); \
  } \
 \
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE A<L, R> get() const { \
    return *this; \
  } \
 \
 private: \
  L _left; \
  R _right; \
}; \
MT_KAHYPAR_DEFINE_VALID_EXPRESSION_BINARY(A); \
 \
template <typename L, typename R, REQUIRE_VALID_WEIGHT_2(L, R)> \
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE auto OPNAME(const L& left, const R& right) { \
  return A<DEDUCE_TYPE(left.get()), DEDUCE_TYPE(right.get())>{left.get(), right.get()}; \
}

// #################  BASIC TYPES  #################

template<typename Subclass>
class HNWeightBase {
 public:
  ~HNWeightBase() {
    ENABLE_IF_DEBUG_WEIGHTS(
      decreaseRefCount();
    )
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Dimension dimension() const {
    return _dimension;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar* get_raw_data() const {
    return _data;
  }

 protected:
  template<typename Other>
  friend class HNWeightBase;

  explicit HNWeightBase(): _data(nullptr), _dimension(0) {}

  explicit HNWeightBase(HNWeightScalar* data, Dimension dimension, uint32_t* ref_count = nullptr):
    _data(data), _dimension(dimension) {
      unused(ref_count);
      ENABLE_IF_DEBUG_WEIGHTS(
        _ref_count = ref_count;
        increaseRefCount();
      )
    }

  HNWeightBase(const HNWeightBase& other):
    _data(other._data),
    _dimension(other._dimension) {
    ENABLE_IF_DEBUG_WEIGHTS(
      _ref_count = other._ref_count;
      increaseRefCount();
    )
  }

  HNWeightBase& operator=(const HNWeightBase& other) {
    assignOther(other);
    return *this;
  }

  HNWeightBase(HNWeightBase&& other):
    _data(other._data),
    _dimension(other._dimension) {
    ENABLE_IF_DEBUG_WEIGHTS(
      _ref_count = other._ref_count;
      other.setEmpty();
    )
  }

  HNWeightBase& operator=(HNWeightBase&& other) {
    assignOther(other);
    return *this;
  }

  template<typename Other>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void assignOther(const HNWeightBase<Other>& other) {
    _data = other._data;
    _dimension = other._dimension;
    ENABLE_IF_DEBUG_WEIGHTS(
      decreaseRefCount();
      _ref_count = other._ref_count;
      increaseRefCount();
    )
  }

  template<typename Other>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void assignOther(HNWeightBase<Other>&& other) {
    _data = other._data;
    _dimension = other._dimension;
    ENABLE_IF_DEBUG_WEIGHTS(
      _ref_count = other._ref_count;
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

 #ifdef MT_KAHYPAR_DEBUG_WEIGHTS
  void setEmpty() {
    decreaseRefCount();
    _data = nullptr;
    _dimension = 0;
    _ref_count = nullptr;
  }

  void increaseRefCount() {
    if (_ref_count != nullptr) {
      __atomic_fetch_add(_ref_count, 1, __ATOMIC_RELAXED);
    }
  }

  void decreaseRefCount() {
    if (_ref_count != nullptr) {
      auto old = __atomic_fetch_sub(_ref_count, 1, __ATOMIC_RELAXED);
      unused(old);
      ALWAYS_ASSERT(old > 0);
    }
  }
 #endif

 protected:
  HNWeightScalar* _data;
  Dimension _dimension;

 #ifdef MT_KAHYPAR_DEBUG_WEIGHTS
  uint32_t* _ref_count = nullptr;
 #endif
};

// ####################  ATOMIC REFERENCE  ####################
class HNWeightAtomicCRef: public HNWeightBase<HNWeightAtomicCRef> {
  using Base = HNWeightBase<HNWeightAtomicCRef>;
  using Base::_data;

 public:
  // construction from references
  template<typename Other>
  explicit HNWeightAtomicCRef(HNWeightBase<Other>&& other, std::memory_order memory_order):
    Base(), _memory_order(memory_order) {
      assignOther(std::move(other));
    }

  // we need copy construction/assignment for atomics so that the get() contract is fulfilled
  HNWeightAtomicCRef(const HNWeightAtomicCRef&) = default;
  HNWeightAtomicCRef& operator=(const HNWeightAtomicCRef&) = default;

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar at(Dimension i) const {
    ASSERT(i < _dimension);
    return __atomic_load_n(&_data[i], transformMemoryOrder(_memory_order));
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightAtomicCRef get() const {
    return HNWeightAtomicCRef(*this, _memory_order);
  }

 private:
  template<typename Other>
  explicit HNWeightAtomicCRef(const HNWeightBase<Other>& other, std::memory_order memory_order):
    Base(other), _memory_order(memory_order) {};

  std::memory_order _memory_order;
};
MT_KAHYPAR_DEFINE_VALID_EXPRESSION(HNWeightAtomicCRef);

// ####################  CONST REFERENCE  ####################
class HNWeightConstRef: public HNWeightBase<HNWeightConstRef> {
  using Base = HNWeightBase<HNWeightConstRef>;

 public:
  explicit HNWeightConstRef(): Base() {}

  explicit HNWeightConstRef(const HNWeightScalar* data, Dimension dimension):
    Base(const_cast<HNWeightScalar*>(data), dimension) {}

  // copy construction/assignment
  HNWeightConstRef(const HNWeightConstRef&) = default;
  HNWeightConstRef& operator=(const HNWeightConstRef&) = default;

  // construction/assignment from non-const reference
  template<typename Other>
  HNWeightConstRef(const HNWeightBase<Other>& other): Base() {
    assignOther(std::move(other));
  }

  template<typename Other>
  HNWeightConstRef& operator=(const HNWeightBase<Other>& other) {
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
class HNWeightRef: public HNWeightBase<HNWeightRef> {
  using Base = HNWeightBase<HNWeightRef>;

  // define structs for atomic expressions
  CREATE_ATOMIC_HNWEIGHTREF_EXPRESSION(FetchAddExpr, __atomic_fetch_add);
  CREATE_ATOMIC_HNWEIGHTREF_EXPRESSION(AddFetchExpr, __atomic_add_fetch);
  CREATE_ATOMIC_HNWEIGHTREF_EXPRESSION(FetchSubExpr, __atomic_fetch_sub);
  CREATE_ATOMIC_HNWEIGHTREF_EXPRESSION(SubFetchExpr, __atomic_sub_fetch);

 public:
  explicit HNWeightRef() : Base() {}

  explicit HNWeightRef(HNWeightScalar* data, Dimension dimension, uint32_t* ref_count = nullptr):
    Base(data, dimension, ref_count) {}

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
      __atomic_store_n(&_data[i], right.at(i), memory_order);
    }
  }

  template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE auto fetch_add(const Expr& right, std::memory_order memory_order) {
    return FetchAddExpr<DEDUCE_TYPE(right.get())>(*this, right.get(), memory_order);
  }

  template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE auto add_fetch(const Expr& right, std::memory_order memory_order) {
    return AddFetchExpr<DEDUCE_TYPE(right.get())>(*this, right.get(), memory_order);
  }

  template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE auto fetch_sub(const Expr& right, std::memory_order memory_order) {
    return FetchSubExpr<DEDUCE_TYPE(right.get())>(*this, right.get(), memory_order);
  }

  template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE auto sub_fetch(const Expr& right, std::memory_order memory_order) {
    return SubFetchExpr<DEDUCE_TYPE(right.get())>(*this, right.get(), memory_order);
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
MT_KAHYPAR_DEFINE_VALID_EXPRESSION_UNARY(HNWeightRef::AddFetchExpr);
MT_KAHYPAR_DEFINE_VALID_EXPRESSION_UNARY(HNWeightRef::FetchSubExpr);
MT_KAHYPAR_DEFINE_VALID_EXPRESSION_UNARY(HNWeightRef::SubFetchExpr);
MT_KAHYPAR_DEFINE_VALID_EXPRESSION(HNWeightRef);

// #################  BROADCAST  #################

class [[nodiscard]] BroadcastExpr {
 public:
  explicit BroadcastExpr(HNWeightScalar value, Dimension dimension): _value(value), _dimension(dimension) { }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar at(Dimension i) const {
    ASSERT(i < _dimension);
    unused(i);
    return _value;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Dimension dimension() const {
    return _dimension;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE BroadcastExpr get() const {
    return *this;
  }

 private:
  HNWeightScalar _value;
  Dimension _dimension;
};
MT_KAHYPAR_DEFINE_VALID_EXPRESSION(BroadcastExpr);

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE BroadcastExpr broadcast(HNWeightScalar value, Dimension dimension) {
  return BroadcastExpr{value, dimension};
}

// ###################  ADDITION & SUBTRACTION  ###################

CREATE_BINARY_HNWEIGHT_OPERATOR(AddExpr, operator+, std::plus<HNWeightScalar>{});
CREATE_BINARY_HNWEIGHT_OPERATOR(SubExpr, operator-, std::minus<HNWeightScalar>{});

// #####################  MIN & MAX  #####################

CREATE_BINARY_HNWEIGHT_OPERATOR(MinExpr, min, std::min);
CREATE_BINARY_HNWEIGHT_OPERATOR(MaxExpr, max, std::max);

// ################  SCALAR MULTIPLICATION & DIVISION  ################

template <typename Scalar, typename R>
class [[nodiscard]] MultExpr {
 public:
  static_assert(is_hypernodeweight_expression<R>);

  explicit MultExpr(Scalar factor, R&& right): _factor(factor), _right(std::move(right)) { }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar at(Dimension i) const {
    return _factor * _right.at(i);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Dimension dimension() const {
    return _right.dimension();
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE MultExpr<Scalar, R> get() const {
    return *this;
  }

 private:
  Scalar _factor;
  R _right;
};

template <typename R, REQUIRE_VALID_WEIGHT(R)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE auto operator*(HNWeightScalar factor, const R& right) {
  return MultExpr<HNWeightScalar, DEDUCE_TYPE(right.get())>{factor, right.get()};
}

template<typename R>
struct IsHypernodeweightExpression<MultExpr<HNWeightScalar, R>> {
  static constexpr bool value = is_hypernodeweight_expression<R>;
};

template <typename R, REQUIRE_VALID_WEIGHT(R)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE auto operator*(double factor, const R& right) {
  return MultExpr<double, DEDUCE_TYPE(right.get())>{factor, right.get()};
}

template<typename R>
struct IsHypernodeweightExpression<MultExpr<double, R>> {
  static constexpr bool value = is_hypernodeweight_expression<R>;
};

template <typename L>
class [[nodiscard]] DivExpr {
 public:
  static_assert(is_hypernodeweight_expression<L>);

  explicit DivExpr(L&& left, HNWeightScalar divisor): _left(std::move(left)), _divisor(divisor) { }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar at(Dimension i) const {
    return _left.at(i) / _divisor;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Dimension dimension() const {
    return _left.dimension();
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE DivExpr<L> get() const {
    return *this;
  }

 private:
  L _left;
  HNWeightScalar _divisor;
};
MT_KAHYPAR_DEFINE_VALID_EXPRESSION_UNARY(DivExpr);

template <typename L, REQUIRE_VALID_WEIGHT(L)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE auto operator/(const L& left, HNWeightScalar divisor) {
  return DivExpr<DEDUCE_TYPE(left.get())>{left, divisor};
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightRef& operator*=(HNWeightRef& weight, HNWeightScalar factor) {
  return weight = factor * weight;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightRef& operator*=(HNWeightRef& weight, double factor) {
  return weight = factor * weight;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightRef& operator/=(HNWeightRef& weight, HNWeightScalar divisor) {
  return weight = weight / divisor;
}

// ###################  COMPARISON  ###################

template <typename L, typename R, REQUIRE_VALID_WEIGHT_2(L, R)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool operator==(const L& left, const R& right) {
  ASSERT(left.dimension() == right.dimension());
  bool is_eq = true;
  for (Dimension i = 0; i < right.dimension(); i++) {
    is_eq &= left.at(i) == right.at(i);
  }
  return is_eq;
}

template <typename L, typename R, REQUIRE_VALID_WEIGHT_2(L, R)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool operator!=(const L& left, const R& right) {
  return !(left == right);
}

template <typename L, typename R, REQUIRE_VALID_WEIGHT_2(L, R)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool operator>(const L& left, const R& right) {
  ASSERT(left.dimension() == right.dimension());
  bool is_greater = true;
  for (Dimension i = 0; i < right.dimension(); i++) {
    is_greater &= left.at(i) > right.at(i);
  }
  return is_greater;
}

template <typename L, typename R, REQUIRE_VALID_WEIGHT_2(L, R)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool operator<=(const L& left, const R& right) {
  ASSERT(left.dimension() == right.dimension());
  bool is_le = true;
  for (Dimension i = 0; i < right.dimension(); i++) {
    is_le &= left.at(i) <= right.at(i);
  }
  return is_le;
}

// TODO: check for bugs related to this
template <typename L, typename R, REQUIRE_VALID_WEIGHT_2(L, R)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool operator<(const L& left, const R& right) {
  ASSERT(left.dimension() == right.dimension());
  bool is_smaller = true;
  for (Dimension i = 0; i < right.dimension(); i++) {
    is_smaller &= left.at(i) < right.at(i);
  }
  return is_smaller;
}

template <typename L, typename R, REQUIRE_VALID_WEIGHT_2(L, R)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool operator>=(const L& left, const R& right) {
  ASSERT(left.dimension() == right.dimension());
  bool is_ge = true;
  for (Dimension i = 0; i < right.dimension(); i++) {
    is_ge &= left.at(i) >= right.at(i);
  }
  return is_ge;
}

// ###################  ARBITRARY MAPPING  ###################

template <typename L, typename Func, bool with_index>
class [[nodiscard]] MappingExpr {
 public:
  static_assert(is_hypernodeweight_expression<L>);
  static_assert(std::is_trivially_copyable_v<Func>);

  explicit MappingExpr(L&& left, Func&& func): _left(std::move(left)), _func(std::forward<Func>(func)) { }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar at(Dimension i) const {
    if constexpr (with_index) {
      return static_cast<HNWeightScalar>(_func(_left.at(i), i));
    } else {
      return static_cast<HNWeightScalar>(_func(_left.at(i)));
    }
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Dimension dimension() const {
    return _left.dimension();
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE MappingExpr<L, Func, with_index> get() const {
    return *this;
  }

 private:
  L _left;
  Func _func;
};

template<typename L, typename Func, bool with_index>
struct IsHypernodeweightExpression<MappingExpr<L, Func, with_index>> {
  static constexpr bool value = is_hypernodeweight_expression<L>;
};

template <typename L, typename Func, REQUIRE_VALID_WEIGHT(L)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE auto map(const L& left, Func&& func) {
  return MappingExpr<DEDUCE_TYPE(left.get()), Func, false>{left.get(), std::forward<Func>(func)};
}

template <typename L, typename Func, REQUIRE_VALID_WEIGHT(L)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE auto mapWithIndex(const L& left, Func&& func) {
  return MappingExpr<DEDUCE_TYPE(left.get()), Func, true>{left.get(), std::forward<Func>(func)};
}

// ##################  UNARY EXPRESSIONS  ##################

template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void eval(Expr&& expr) {
  for (Dimension i = 0; i < expr.dimension(); i++) {
    expr.at(i);
  }
}

template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar sum(const Expr& expr) {
  HNWeightScalar sum = 0;
  for (Dimension i = 0; i < expr.dimension(); i++) {
    sum += expr.at(i);
  }
  return sum;
}

template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar squaresSum(const Expr& expr) {
  HNWeightScalar sum = 0;
  for (Dimension i = 0; i < expr.dimension(); i++) {
    sum += expr.at(i) * expr.at(i);
  }
  return sum;
}

template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar maxOf(const Expr& expr) {
  HNWeightScalar max = std::numeric_limits<HNWeightScalar>::min();
  for (Dimension i = 0; i < expr.dimension(); i++) {
    max = std::max(max, expr.at(i));
  }
  return max;
}

template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar minOf(const Expr& expr) {
  HNWeightScalar min = std::numeric_limits<HNWeightScalar>::max();
  for (Dimension i = 0; i < expr.dimension(); i++) {
    min = std::min(min, expr.at(i));
  }
  return min;
}

template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool isZero(const Expr& expr) {
  bool is_zero = true;
  for (Dimension i = 0; i < expr.dimension(); i++) {
    is_zero &= expr.at(i) == 0;
  }
  return is_zero;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightConstRef newInvalid() {
  return HNWeightConstRef();
}

template <typename Other>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool isInvalid(const HNWeightBase<Other>& value) {
  return value.dimension() == 0;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightConstRef toNonAtomic(const HNWeightAtomicCRef& value) {
  return HNWeightConstRef(value.get_raw_data(), value.dimension());
}

template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE std::string toString(const Expr& expr) {
  std::stringstream ss;
  ss << expr;
  return ss.str();
}

// ###################  PRINT  ###################

template <typename R, REQUIRE_VALID_WEIGHT(R)>
std::ostream & operator<< (std::ostream& os, const R& right) {
  if (right.dimension() == 1) {
    return os << right.at(0);
  }

  os << "[ ";
  for (Dimension i = 0; i < right.dimension(); i++) {
    os << std::to_string(right.at(i)) + " ";
  }
  return os << "]";
}

// ###################  ALLOCATE  ###################

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

} // namespace weight

using weight::Dimension;
using weight::HNWeightScalar;
using weight::HNWeightRef;
using weight::HNWeightConstRef;
using weight::HNWeightAtomicCRef;

} // namespace mt_kahypar

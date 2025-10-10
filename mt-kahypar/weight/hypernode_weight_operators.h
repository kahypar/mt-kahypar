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

#include <limits>
#include <ostream>
#include <sstream>

#include "mt-kahypar/macros.h"
#include "mt-kahypar/weight/hypernode_weight_base.h"

namespace mt_kahypar {
namespace weight {

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
  return DivExpr<DEDUCE_TYPE(left.get())>{left.get(), divisor};
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightRef& operator*=(HNWeightRef& weight, HNWeightScalar factor) {
  return weight = factor * weight;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightRef& operator*=(HNWeightRef&& weight, HNWeightScalar factor) {
  return weight = factor * weight;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightRef& operator*=(HNWeightRef& weight, double factor) {
  return weight = factor * weight;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightRef& operator*=(HNWeightRef&& weight, double factor) {
  return weight = factor * weight;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightRef& operator/=(HNWeightRef& weight, HNWeightScalar divisor) {
  return weight = weight / divisor;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightRef& operator/=(HNWeightRef&& weight, HNWeightScalar divisor) {
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


// ###################  TERNARY OPERATOR  ###################

template <typename L, typename R>
class [[nodiscard]] TernaryExpr {
 public:
  static_assert(is_hypernodeweight_expression<L> && is_hypernodeweight_expression<R>);

  explicit TernaryExpr(bool condition, L&& left, R&& right): _condition(condition), _left(std::move(left)), _right(std::move(right)) {
    ASSERT(left.dimension() == right.dimension());
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar at(Dimension i) const {
    return _condition ? _left.at(i) : _right.at(i);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Dimension dimension() const {
    return _left.dimension();
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE TernaryExpr<L, R> get() const {
    return *this;
  }

 private:
  bool _condition;
  L _left;
  R _right;
};
MT_KAHYPAR_DEFINE_VALID_EXPRESSION_BINARY(TernaryExpr);

template <typename L, typename R, REQUIRE_VALID_WEIGHT_2(L, R)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE auto ternary(bool condition, const L& left, const R& right) {
  return TernaryExpr<DEDUCE_TYPE(left.get()), DEDUCE_TYPE(right.get())>{condition, left.get(), right.get()};
}

template <typename Func>
class [[nodiscard]] LazyExpr {
  static_assert(std::is_trivially_copyable_v<Func>);

 public:
  explicit LazyExpr(Func&& func, Dimension dimension):
    _func(std::forward<Func>(func)), _dimension(dimension) { }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightScalar at(Dimension i) const {
    return _func().at(i);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Dimension dimension() const {
    return _dimension;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE LazyExpr<Func> get() const {
    return *this;
  }

 private:
  Func _func;
  Dimension _dimension;
};

template<typename Func>
struct IsHypernodeweightExpression<LazyExpr<Func>> {
  static constexpr bool value = true;
};

template <typename Func>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE LazyExpr<Func> lazy(Func&& func, Dimension dimension) {
  return LazyExpr<Func>{std::forward<Func>(func), dimension};
}


// ##################  UNARY EXPRESSIONS  ##################

template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE auto operator-(const Expr& expr) {
  return map(expr, [](HNWeightScalar val) { return -val; });
}

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

template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool isInvalid(const Expr& expr) {
  return expr.dimension() == 0;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HNWeightConstRef toNonAtomic(const HNWeightAtomicCRef& value) {
  return HNWeightConstRef(value.get_raw_data(), value.dimension());
}

template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Expr copy(const Expr& value) {
  return Expr(value.get_raw_data(), value.dimension());
}

template <typename Other>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void replace(HNWeightRefBase<Other>& lhs, const HNWeightRefBase<Other>& value) {
  lhs.set_raw_data(value.get_raw_data(), value.dimension());
}


// ###################  PRINTING  ###################

template <typename R, REQUIRE_VALID_WEIGHT(R)>
std::ostream& operator<<(std::ostream& os, const R& right) {
  return os << toString(right);
}

template <typename Expr, REQUIRE_VALID_WEIGHT(Expr)>
std::string toString(const Expr& expr) {
  if (expr.dimension() == 1) {
    return std::to_string(expr.at(0));
  }

  std::stringstream ss;
  ss << "[ ";
  for (Dimension i = 0; i < expr.dimension(); i++) {
    ss << expr.at(i) << " ";
  }
  ss << "]";
  return ss.str();
}

} // namespace weight
} // namespace mt_kahypar

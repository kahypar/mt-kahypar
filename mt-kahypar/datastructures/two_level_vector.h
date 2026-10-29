#pragma once

#include <vector>
#include <unordered_map>
#include <limits>
#include <algorithm>

namespace mt_kahypar { namespace ds {

// Two-level vector: compact base array (uint8_t/uint16_t) + sparse overflow map
// Default value semantics: if base is empty, get(i) returns 1.
template <typename BaseT, typename ValueT>
class TwoLevelVector {
 public:
  TwoLevelVector() = default;

  bool empty() const { return _base.empty(); }
  size_t size() const { return _base.size(); }

  void init(size_t n, ValueT default_value) {
    _base.assign(n, clamp_to_base(default_value));
    _overflow.clear();
    if (default_value > max_base()) {
      const ValueT rem = default_value - static_cast<ValueT>(max_base());
      for (size_t i = 0; i < n; ++i) _overflow.emplace(i, rem);
    }
  }

  // Lazy init with default = 1
  void ensure_initialized(size_t n) {
    if (_base.empty()) _base.assign(n, clamp_to_base(1));
  }

  ValueT get(size_t i) const {
    if (_base.empty()) return static_cast<ValueT>(1);
    ValueT lo = static_cast<ValueT>(_base[i]);
    auto it = _overflow.find(i);
    if (it != _overflow.end()) lo = static_cast<ValueT>(lo + it->second);
    return lo;
  }

  void set(size_t i, ValueT v) {
    if (v == static_cast<ValueT>(1) && _base.empty()) {
      // Keep lazy default without allocating
      return;
    }
    if (_base.empty()) ensure_initialized(i + 1);
    if (_base.size() <= i) _base.resize(i + 1, clamp_to_base(1));
    if (v <= max_base()) {
      _base[i] = static_cast<BaseT>(v);
      // erase any overflow
      auto it = _overflow.find(i);
      if (it != _overflow.end()) _overflow.erase(it);
    } else {
      _base[i] = max_base();
      _overflow[i] = static_cast<ValueT>(v - max_base());
    }
  }

  // Writable reference proxy returned by non-const operator[]
  class Ref {
   public:
    Ref(TwoLevelVector& owner, size_t idx) : _owner(owner), _idx(idx) {}
    // Implicit read
    operator ValueT() const { return _owner.get(_idx); }
    // Assignment write
    Ref& operator=(const ValueT v) {
      _owner.set(_idx, v);
      return *this;
    }
    // Support assigning from another Ref
    Ref& operator=(const Ref& other) {
      return (*this = static_cast<ValueT>(other));
    }
   private:
    TwoLevelVector& _owner;
    size_t _idx;
  };

  // Read-only access
  ValueT operator[](size_t i) const { return get(i); }
  // Read-write access via proxy
  Ref operator[](size_t i) { return Ref(*this, i); }

  void clear() {
    _base.clear();
    _base.shrink_to_fit();
    _overflow.clear();
  }

  // Approximate memory usage in bytes (does not include hash table overhead)
  size_t approx_bytes() const {
    return sizeof(BaseT) * _base.capacity() + _overflow.size() * (sizeof(size_t) + sizeof(ValueT));
  }

 private:
  static constexpr BaseT max_base() { return std::numeric_limits<BaseT>::max(); }
  static BaseT clamp_to_base(ValueT v) {
    return static_cast<BaseT>(std::min<ValueT>(v, max_base()));
  }

  std::vector<BaseT> _base;                         // compact base
  std::unordered_map<size_t, ValueT> _overflow;     // sparse remainder per index
};

}} // namespace mt_kahypar::ds

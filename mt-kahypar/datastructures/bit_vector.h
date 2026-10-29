/*******************************************************************************
 * Dynamic Bit Vector
 * A minimal, dependency-free dynamic bitset used in Mt-KaHyPar.
 * Stores bits in a std::vector<uint64_t> and provides vector<bool>-like access.
 ******************************************************************************/

#pragma once

#include <vector>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <stdexcept>

namespace mt_kahypar {
namespace ds {

class BitVector {
 public:
  BitVector() : _size_bits(0), _words() {}

  // Copy/move
  BitVector(const BitVector&) = default;
  BitVector(BitVector&&) noexcept = default;
  BitVector& operator=(const BitVector&) = default;
  BitVector& operator=(BitVector&&) noexcept = default;

  // A proxy reference similar to std::vector<bool>::reference
  class reference {
   public:
    reference(uint64_t& word, uint64_t mask) : _word(word), _mask(mask) {}
    reference& operator=(bool value) {
      if (value) { _word |= _mask; } else { _word &= ~_mask; }
      return *this;
    }
    reference& operator=(const reference& other) {
      return (*this = static_cast<bool>(other));
    }
    operator bool() const { return (_word & _mask) != 0; }
    void flip() { _word ^= _mask; }
   private:
    uint64_t& _word;
    const uint64_t _mask;
  };

  // Const element access returns a bool
  bool operator[](size_t idx) const {
    if (idx >= _size_bits) {
      throw std::out_of_range("BitVector index out of range");
    }
    const auto wi = word_index(idx);
    const uint64_t mask = bit_mask(idx);
    return (wi < _words.size()) && ((_words[wi] & mask) != 0);
  }

  // Mutable element access returns proxy
  reference operator[](size_t idx) {
    ensure_index(idx);
    const auto wi = word_index(idx);
    const uint64_t mask = bit_mask(idx);
    return reference(_words[wi], mask);
  }

  // Append one bit
  void push_back(bool value) {
    const size_t idx = _size_bits;
    // Ensure underlying storage has capacity for this bit but do not
    // change logical size until after writing the bit.
    const size_t wi = word_index(idx);
    if (wi >= _words.size()) {
      _words.resize(wi + 1, 0ULL);
    }
    if (value) {
      _words[wi] |= bit_mask(idx);
    }
    ++_size_bits;
  }

  // Assign n bits all to value
  void assign(size_t n, bool value) {
    _size_bits = n;
    const size_t needed_words = words_for_bits(_size_bits);
    _words.assign(needed_words, value ? ~uint64_t(0) : uint64_t(0));
    // Mask off unused high bits in last word
    if (needed_words > 0) {
      const size_t rem = (_size_bits & 63ULL);
      if (rem != 0) {
        const uint64_t mask = (uint64_t(1) << rem) - 1ULL;
        _words.back() &= mask;
      }
    }
  }

  void clear() {
    _size_bits = 0;
    _words.clear();
  }

  void shrink_to_fit() { _words.shrink_to_fit(); }

  bool empty() const { return _size_bits == 0; }
  size_t size() const { return _size_bits; }

  // Capacity in bits (based on underlying vector capacity)
  size_t capacity() const { return _words.capacity() * 64ULL; }

  // Reserve capacity for at least nbits bits without changing size
  void reserve(size_t nbits) { _words.reserve(words_for_bits(nbits)); }

  // Raw word storage access (const)
  const std::vector<uint64_t>& words() const { return _words; }

 private:
  static size_t words_for_bits(size_t nbits) { return (nbits + 63ULL) / 64ULL; }
  static size_t word_index(size_t bit_index) { return bit_index / 64ULL; }
  static uint64_t bit_mask(size_t bit_index) { return 1ULL << (bit_index & 63ULL); }

  void ensure_index(size_t bit_index) {
    const size_t wi = word_index(bit_index);
    if (wi >= _words.size()) {
      _words.resize(wi + 1, 0ULL);
    }
    if (bit_index >= _size_bits) {
      _size_bits = bit_index + 1; // logical size grows on write access
    }
  }

  size_t _size_bits;
  std::vector<uint64_t> _words;
};

} // namespace ds
} // namespace mt_kahypar

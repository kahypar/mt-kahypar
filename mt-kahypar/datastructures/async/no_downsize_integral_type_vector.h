//
// Created by mlaupichler on 31.08.21.
//

#pragma once

namespace mt_kahypar::ds {

    // ! Vector based data structure for primitive integral types that does not free memory until destruction in order to avoid re-allocations as much
    // ! as possible. Provides amortized O(1) operations at back and front as well as O(1) clear.
    template<typename T>
    class NoDownsizeIntegralTypeVector {

        static_assert(std::is_integral<T>::value,
                      "Type T for NoDownsizeIntegralTypeVector has to be an integral type!");

    public:

        NoDownsizeIntegralTypeVector() :
            _num_expired_threshold(1),
            _first_valid(0),
            _first_invalid(0),
            _data() {}

        explicit NoDownsizeIntegralTypeVector(const size_t initial_capacity, const T& initial_value = std::numeric_limits<T>::max()) :
            _num_expired_threshold(std::max(initial_capacity / 4, size_t(1))),
            _first_valid(0),
            _first_invalid(0),
            _data(initial_capacity, initial_value) {}

        NoDownsizeIntegralTypeVector(const NoDownsizeIntegralTypeVector& other) = default;
        NoDownsizeIntegralTypeVector& operator=(const NoDownsizeIntegralTypeVector& other) = default;
        NoDownsizeIntegralTypeVector(NoDownsizeIntegralTypeVector&& other)  noexcept = default;
        NoDownsizeIntegralTypeVector& operator=(NoDownsizeIntegralTypeVector&& other)  noexcept = default;

        size_t size() const {
          return _first_invalid - _first_valid;
        }

        bool empty() const {
          return _first_invalid == _first_valid;
        }

        size_t size_in_bytes() const {
          return (size_t) _data.size() * sizeof(T);
        };

        void push_back(const T element) {
          if (_first_invalid == _data.size()) {
            if (_first_valid >= _num_expired_threshold) {
              // If at least as many elements as the threshold are expired, remove expired part
              remove_expired();
            } else {
              // Else have std::vector deal with the push_back on full vector and return
              _data.push_back(element);
              _num_expired_threshold = _data.size() / 4;
              ++_first_invalid;
              return;
            }
          }
          ASSERT(_first_invalid < _data.size());
          _data[_first_invalid] = element;
          ++_first_invalid;
        }

        T pop_back() {
          ASSERT(_first_invalid > _first_valid);
          const T back = _data[--_first_invalid];
          if (_first_invalid == _first_valid) {
            // Remove expired in O(1) if valid edges became empty
            remove_expired();
          }
          return back;
        }

        T pop_front() {
          ASSERT(_first_invalid > _first_valid);
          const T front = _data[_first_valid++];
          if (_first_invalid == _first_valid) {
            // Remove expired in O(1) if valid edges became empty
            remove_expired();
          }
          return front;
        }

        T& front() {
          ASSERT(_first_invalid > _first_valid);
          return _data[_first_valid];
        }

        T& back() {
          ASSERT(_first_invalid > _first_valid);
          return _data[_first_invalid - 1];
        }

        // ! Clear in O(1). Has no effect on capacity, i.e. does not free any memory previously allocated on construction or push_back()
        void clear() {
          ASSERT(_first_invalid >= _first_valid);
          _first_valid = 0;
          _first_invalid = 0;
        }

        auto begin() {ASSERT(_first_valid <= _data.size()); return _data.begin() + _first_valid;};
        auto end() {ASSERT(_first_invalid <= _data.size()); return _data.begin() + _first_invalid;};
        auto cbegin() const {ASSERT(_first_valid <= _data.size()); return _data.cbegin() + _first_valid;};
        auto cend() const {ASSERT(_first_invalid <= _data.size()); return _data.cbegin() + _first_invalid;};
        auto begin() const {ASSERT(_first_valid <= _data.size()); return _data.begin() + _first_valid;};
        auto end() const {ASSERT(_first_invalid <= _data.size()); return _data.begin() + _first_invalid;};

    private:

        void remove_expired() {
          for (size_t i = 0; i < size(); ++i) {
            _data[i] = _data[i + _first_valid];
          }
          _first_invalid -= _first_valid;
          _first_valid = 0;
        }

        size_t _num_expired_threshold;
        // Marks first index of the currently valid elements (everything before this is considered expired)
        size_t _first_valid;
        // Marks first index past the currently valid elements
        size_t _first_invalid;
        vec<T> _data;
    };

}  // end namespace

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

#include <numeric>
#include <vector>

#include <tbb/parallel_scan.h>

#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {

  template<typename InIt, typename OutIt, typename BinOp>
  struct ParallelPrefixSumBody {
    using T = typename ::std::iterator_traits<InIt>::value_type;

    InIt first;
    OutIt out;
    T sum, neutral_element;
    BinOp& f;

    ParallelPrefixSumBody(InIt first, OutIt out, T neutral_element, BinOp& f):
      first(first),
      out(out),
      sum(neutral_element),
      neutral_element(neutral_element),
      f(f) { }

    ParallelPrefixSumBody(ParallelPrefixSumBody& other, tbb::split) :
      first(other.first),
      out(other.out),
      sum(other.neutral_element),
      neutral_element(other.neutral_element),
      f(other.f) { }

    void operator()(const tbb::blocked_range<size_t>& r, tbb::pre_scan_tag ) {
      for (size_t i = r.begin(); i < r.end(); ++i) {
        sum = f(sum, *(first + i));
      }
    }

    void operator()(const tbb::blocked_range<size_t>& r, tbb::final_scan_tag ) {
      for (size_t i = r.begin(); i < r.end(); ++i) {
        sum = f(sum, *(first + i));
        *(out + i) = sum;
      }
    }

    void reverse_join(ParallelPrefixSumBody& other) {
      sum = f(sum, other.sum);
    }

    void assign(ParallelPrefixSumBody& other) {
      sum = other.sum;
    }

  };

  template <class InIt, class OutIt, class BinOp>
  static void sequential_prefix_sum(InIt first, InIt last, OutIt d, typename std::iterator_traits<InIt>::value_type init, BinOp f) {
    while (first != last) {
      init = f(init, *first);
      *d = init;
      ++d;
      ++first;
    }
  }

  template <class InIt, class OutIt, class BinOp>
  static void parallel_prefix_sum(InIt first, InIt last, OutIt d, BinOp f,
                                  typename std::iterator_traits<InIt>::value_type neutral_element) {

    typename std::iterator_traits<InIt>::difference_type n = last - first;

    if (n < (1 << 16)) {
      return sequential_prefix_sum(first, last, d, neutral_element, f);
    }

    ParallelPrefixSumBody<InIt, OutIt, BinOp> body(first, d, neutral_element, f);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0, static_cast<size_t>(n)), body);
  }

}

namespace mt_kahypar::parallel {

template<typename T,
         template<class> class V = parallel::scalable_vector>
class TBBPrefixSum {

 public:
  TBBPrefixSum(V<T>& data) :
    _sum(0),
    _data(data) { }

  TBBPrefixSum(TBBPrefixSum& prefix_sum, tbb::split) :
    _sum(0),
    _data(prefix_sum._data) { }

  T total_sum() const {
    return _sum;
  }

  size_t size() const {
    return _data.size() + 1;
  }

  T operator[] (const size_t i) const {
    ASSERT(i <= _data.size());
    if ( i > 0 ) {
      return _data[i - 1];
    } else {
      return static_cast<T>(0);
    }
  }

  T value(const size_t i) const {
    ASSERT(i < _data.size(), V(i) << V(_data.size()));
    if ( i > 0 ) {
      return _data[i] - _data[i - 1];
    } else {
      return _data[0];
    }
  }

  template<typename Tag>
  void operator()(const tbb::blocked_range<size_t>& range, Tag) {
      T temp = _sum;
      for( size_t i = range.begin(); i < range.end(); ++i ) {
          temp = temp + _data[i];
          if( Tag::is_final_scan() ) {
            _data[i] = temp;
          }
      }
      _sum = temp;
  }

  void reverse_join(TBBPrefixSum& prefix_sum) {
    _sum += prefix_sum._sum;
  }

  void assign(TBBPrefixSum& prefix_sum) {
    _sum = prefix_sum._sum;
  }

 private:
  T _sum;
  V<T>& _data;
};

}  // namespace mt_kahypar::parallel

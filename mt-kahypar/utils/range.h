/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include <vector>

template<typename IteratorT>
class IteratorRange {
public:
  IteratorRange(const IteratorT& first, const IteratorT& firstInvalid) : __begin(first), __end(firstInvalid) { }

  using Iterator = IteratorT; // make publicly visible

  IteratorT begin() {
    return __begin;
  }

  IteratorT end() {
    return __end;
  }

  bool empty() {
    return __begin == __end;
  }

private:
  IteratorT __begin, __end;
};



template<typename RangeT>
class ConcatenatedRange {
private:
  struct begin_tag {};
  struct end_tag {};

public:
  class Iterator {
  public:
    Iterator(std::vector<RangeT>& ranges, begin_tag) : ranges(ranges), currentRange(0), currentRangeIterator(ranges.front().begin()) {
      moveToNextRange();
    }

    Iterator(std::vector<RangeT>& ranges, end_tag) : ranges(ranges), currentRange(ranges.size() - 1), currentRangeIterator(ranges.back().end()) { }


    bool operator==(Iterator& o) {
      return currentRange == o.currentRange && currentRangeIterator == o.currentRangeIterator;
    }

    bool operator!=(Iterator& o) {
      return !operator==(o);
    }

    Iterator& operator++() {
      // if we're at the end of the current range, advance to the next
      // restrict currentRange to ranges.size() - 1, since the end() ConcatenatedRange::Iterator has to be initialized somehow
      if (++currentRangeIterator == ranges[currentRange].end()) {
        moveToNextRange();
      }
      return *this;
    }

    typename RangeT::Iterator::value_type operator*() const {
      return *currentRangeIterator;
    }

  private:
    std::vector<RangeT>& ranges;
    size_t currentRange;
    typename RangeT::Iterator currentRangeIterator;

    void moveToNextRange() {
      while (currentRangeIterator == ranges[currentRange].end() && currentRange < ranges.size() - 1) {
        currentRangeIterator = ranges[++currentRange].begin();
      }
    }

  };

  Iterator begin() {
    assert(!ranges.empty());
    return Iterator(ranges, begin_tag());
  }

  Iterator end() {
    assert(!ranges.empty());
    return Iterator(ranges, end_tag());
  }

  void concat(RangeT&& r) {
    ranges.push_back(r);
  }

  void concat(RangeT& r) {
    ranges.push_back(r);
  }

private:
  std::vector<RangeT> ranges;
};

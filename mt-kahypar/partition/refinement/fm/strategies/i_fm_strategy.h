/*******************************************************************************
 * MIT License
 *
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include "mt-kahypar/macros.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/cast.h"

namespace mt_kahypar {

// TODO: this is still a bit hacky, is there any better way?
struct localized_k_way_fm_s;

struct localized_k_way_fm_t {
  localized_k_way_fm_s* local_fm;
  mt_kahypar_partition_type_t type;
};

namespace utils {
// compare cast.h
template<typename LocalFM>
localized_k_way_fm_t localized_fm_cast(LocalFM& local_fm) {
  return localized_k_way_fm_t {
    reinterpret_cast<localized_k_way_fm_s*>(&local_fm), LocalFM::PartitionedHypergraph::TYPE };
}

template<typename LocalFM>
LocalFM& cast(localized_k_way_fm_t fm) {
  if ( LocalFM::PartitionedHypergraph::TYPE != fm.type ) {
    ERR("Cannot cast local FM [" << typeToString(fm.type) << "to"
        << typeToString(LocalFM::PartitionedHypergraph::TYPE) << "]");
  }
  return *reinterpret_cast<LocalFM*>(fm.local_fm);
}

} // namespace utils


class IFMStrategy {
 public:
  IFMStrategy(const IFMStrategy&) = delete;
  IFMStrategy(IFMStrategy&&) = delete;
  IFMStrategy & operator= (const IFMStrategy &) = delete;
  IFMStrategy & operator= (IFMStrategy &&) = delete;

  virtual ~IFMStrategy() = default;

  bool dispatchedFindMoves(localized_k_way_fm_t local_fm, mt_kahypar_partitioned_hypergraph_t& phg,
                           size_t task_id, size_t num_seeds, size_t round) {
    return dispatchedFindMovesImpl(local_fm, phg, task_id, num_seeds, round);
  }

  bool isUnconstrainedRound(size_t round) const {
    return isUnconstrainedRoundImpl(round);
  }

  bool includesUnconstrained() const {
    return includesUnconstrainedImpl();
  }

  void reportImprovement(size_t round, Gain absolute_improvement, double relative_improvement) {
    reportImprovementImpl(round, absolute_improvement, relative_improvement);
  }

 protected:
  IFMStrategy() = default;

 private:
  virtual bool dispatchedFindMovesImpl(localized_k_way_fm_t local_fm, mt_kahypar_partitioned_hypergraph_t& phg,
                                       size_t task_id, size_t num_seeds, size_t round) = 0;

  virtual bool isUnconstrainedRoundImpl(size_t round) const = 0;

  virtual bool includesUnconstrainedImpl() const = 0;

  virtual void reportImprovementImpl(size_t, Gain, double) {
    // most strategies don't use this
  }
};

}  // namespace mt_kahypar

/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/macros.h"

namespace mt_kahypar {

class CutRollback {

 public:
  struct RecalculationData {
    MoveID first_out, last_in;
    HypernodeID moved_out;
    RecalculationData() :
      first_out(std::numeric_limits<MoveID>::max()),
      last_in(std::numeric_limits<MoveID>::min()),
      moved_out(0)
      { }

    void reset() {
      first_out = std::numeric_limits<MoveID>::max();
      last_in = std::numeric_limits<MoveID>::min();
      moved_out = 0;
    }
  };

  static void updateMove(const MoveID m_id,
                         const Move& m,
                         vec<RecalculationData>& r) {
    r[m.from].first_out = std::min(r[m.from].first_out, m_id);
    r[m.to].last_in = std::max(r[m.to].last_in, m_id);
    ++r[m.from].moved_out;
  }

  static void updateNonMovedPinInBlock(const PartitionID,
                                       vec<RecalculationData>&) {
    // Do nothing here
  }

  template<typename PartitionedHypergraph>
  static bool hasBenefit(const PartitionedHypergraph& phg,
                         const HyperedgeID e,
                         const MoveID m_id,
                         const Move& m,
                         vec<RecalculationData>& r) {
    const HypernodeID edge_size = phg.edgeSize(e);
    const bool was_potentially_non_cut_edge_at_some_point =
      phg.pinCountInPart(e, m.to) + r[m.to].moved_out == edge_size;
    return was_potentially_non_cut_edge_at_some_point && r[m.to].last_in == m_id && m_id < r[m.to].first_out;
  }

  template<typename PartitionedHypergraph>
  static bool hasPenalty(const PartitionedHypergraph& phg,
                         const HyperedgeID e,
                         const MoveID m_id,
                         const Move& m,
                         vec<RecalculationData>& r) {
    const HypernodeID edge_size = phg.edgeSize(e);
    const bool was_potentially_non_cut_edge_at_some_point =
      phg.pinCountInPart(e, m.from) + r[m.from].moved_out == edge_size;
    return was_potentially_non_cut_edge_at_some_point && r[m.from].first_out == m_id && m_id > r[m.from].last_in;
  }
};

}  // namespace mt_kahypar

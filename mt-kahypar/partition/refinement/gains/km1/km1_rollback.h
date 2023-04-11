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

#include <algorithm>

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/macros.h"

namespace mt_kahypar {

class Km1Rollback {

 public:
  struct RecalculationData {
    MoveID first_in, last_out;
    HypernodeID remaining_pins;
    RecalculationData() :
      first_in(std::numeric_limits<MoveID>::max()),
      last_out(std::numeric_limits<MoveID>::min()),
      remaining_pins(0)
      { }

    void reset() {
      first_in = std::numeric_limits<MoveID>::max();
      last_out = std::numeric_limits<MoveID>::min();
      remaining_pins = 0;
    }
  };

  static void updateMove(const MoveID m_id,
                         const Move& m,
                         vec<RecalculationData>& r) {
    r[m.to].first_in = std::min(r[m.to].first_in, m_id);
    r[m.from].last_out = std::max(r[m.from].last_out, m_id);
  }

  static void updateNonMovedPinInBlock(const PartitionID block,
                                       vec<RecalculationData>& r) {
    r[block].remaining_pins++;
  }

  template<typename PartitionedHypergraph>
  static bool hasBenefit(const PartitionedHypergraph&,
                         const HyperedgeID,
                         const MoveID m_id,
                         const Move& m,
                         vec<RecalculationData>& r) {
    return r[m.from].last_out == m_id && r[m.from].first_in > m_id && r[m.from].remaining_pins == 0;
  }

  template<typename PartitionedHypergraph>
  static bool hasPenalty(const PartitionedHypergraph&,
                         const HyperedgeID,
                         const MoveID m_id,
                         const Move& m,
                         vec<RecalculationData>& r) {
    return r[m.to].first_in == m_id && r[m.to].last_out < m_id && r[m.to].remaining_pins == 0;
  }
};

}  // namespace mt_kahypar

/*******************************************************************************
 * This file is part of MT-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include <mt-kahypar/datastructures/hypergraph_common.h>

namespace mt_kahypar {

// adaptive random walk stopping rule from KaHyPar
class StopRule {
public:
  StopRule(HypernodeID numNodes) : beta(std::log(numNodes)) { }

  bool searchShouldStop() {
    return (numSteps > beta) && (Mk == 0 || numSteps >= ( variance / (Mk*Mk) ) * stopFactor );
  }

  void update(Gain gain) {
    ++numSteps;
    if (numSteps == 1) {
      Mk = gain;
      MkPrevious = gain;
      SkPrevious = 0.0;
      variance = 0.0;
    } else {
      Mk = MkPrevious + (gain - MkPrevious) / numSteps;
      Sk = SkPrevious + (gain - MkPrevious) * (gain - Mk);
      variance = Sk / (numSteps - 1.0);

      MkPrevious = Mk;
      SkPrevious = Sk;
    }
  }

  void reset() {
    numSteps = 0;
    variance = 0.0;
  }

private:
  size_t numSteps = 0;
  double variance = 0.0, Mk = 0.0, MkPrevious = 0.0, Sk = 0.0, SkPrevious = 0.0;
  const double alpha = 1.0;   // make parameter if it doesn't work well
  const double stopFactor = (alpha / 2.0) - 0.25;
  double beta;
};
}
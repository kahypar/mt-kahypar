/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {
class LastRatingWins {
 public:
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static bool acceptEqual(const int) {
    return true;
  }

  LastRatingWins(const LastRatingWins&) = delete;
  LastRatingWins & operator= (const LastRatingWins &) = delete;

  LastRatingWins(LastRatingWins&&) = delete;
  LastRatingWins & operator= (LastRatingWins &&) = delete;

 protected:
  ~LastRatingWins() = default;
};

class FirstRatingWins {
 public:
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static bool acceptEqual(const int) {
    return false;
  }

  FirstRatingWins(const FirstRatingWins&) = delete;
  FirstRatingWins & operator= (const FirstRatingWins &) = delete;

  FirstRatingWins(FirstRatingWins&&) = delete;
  FirstRatingWins & operator= (FirstRatingWins &&) = delete;

 protected:
  ~FirstRatingWins() = default;
};

class RandomRatingWins {
 public:
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static bool acceptEqual(const int cpu_id) {
    return utils::Randomize::instance().flipCoin(cpu_id);
  }

  RandomRatingWins(const RandomRatingWins&) = delete;
  RandomRatingWins & operator= (const RandomRatingWins &) = delete;

  RandomRatingWins(RandomRatingWins&&) = delete;
  RandomRatingWins & operator= (RandomRatingWins &&) = delete;

 protected:
  ~RandomRatingWins() = default;
};
}  // namespace mt_kahypar

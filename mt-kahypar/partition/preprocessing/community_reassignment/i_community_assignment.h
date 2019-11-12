/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/definitions.h"

namespace mt_kahypar {
namespace preprocessing {

class ICommunityAssignment {
 public:
  ICommunityAssignment(const ICommunityAssignment&) = delete;
  ICommunityAssignment& operator= (const ICommunityAssignment&) = delete;
  ICommunityAssignment(ICommunityAssignment&&) = delete;
  ICommunityAssignment& operator= (ICommunityAssignment&&) = delete;

  virtual ~ICommunityAssignment() = default;

  std::vector<PartitionID> computeAssignment() {
    return computeAssignmentImpl();
  }

 protected:
  ICommunityAssignment() = default;

 private:
  virtual std::vector<PartitionID> computeAssignmentImpl() = 0;
};

} // namespace preprocessing
} // namespace mt_kahypar
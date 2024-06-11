/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2017 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2017 Robin Andre <robinandre1995@web.de>
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
#include "individual.h"

namespace mt_kahypar {
std::ostream& operator<< (std::ostream& os, const Individual& individual) {
  os << "Fitness: " << individual.fitness() << std::endl;
  os << "Partition:------------------------------------" << std::endl;
  for (size_t i = 0; i < individual.partition().size(); ++i) {
    os << individual.partition()[i] << " ";
  }
  return os;
}
}  // namespace mt_kahypar
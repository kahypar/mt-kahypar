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

#include <iostream>
#include <string>
  
namespace mt_kahypar {

enum class Type : int8_t {
  Unweighted = 0,
  EdgeWeights = 1,
  NodeWeights = 10,
  EdgeAndNodeWeights = 11,
};

std::ostream& operator<< (std::ostream& os, const Type& type) {
  switch (type) {
    case Type::Unweighted: return os << "unweighted";
    case Type::EdgeWeights: return os << "edge_weights";
    case Type::NodeWeights: return os << "node_weights";
    case Type::EdgeAndNodeWeights: return os << "edge_and_node_weights";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(type);
}


} // namesapce mt_kahypar
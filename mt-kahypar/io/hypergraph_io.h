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

#include <cstring>
#include <fstream>
#include <iostream>
#include <thread>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>

#include "tbb/parallel_for.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
namespace io {
  Hypergraph readHypergraphFile(const std::string& filename,
                                const TaskGroupID task_group_id,
                                const bool stable_construction_of_incident_edges = false);
  void readPartitionFile(const std::string& filename, std::vector<PartitionID>& partition);
  void writePartitionFile(const PartitionedHypergraph& phg, const std::string& filename);

}  // namespace io
}  // namespace mt_kahypar

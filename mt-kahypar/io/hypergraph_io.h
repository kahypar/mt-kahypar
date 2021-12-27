/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <string>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace io {

  using Hyperedge = parallel::scalable_vector<HypernodeID>;
  using HyperedgeVector = parallel::scalable_vector<Hyperedge>;
  #ifdef USE_GRAPH_PARTITIONER
  using EdgeVector = parallel::scalable_vector<std::pair<HypernodeID, HypernodeID>>;
  #else
  using EdgeVector = HyperedgeVector;
  #endif

  void readHypergraphFile(const std::string& filename,
                          HyperedgeID& num_hyperedges,
                          HypernodeID& num_hypernodes,
                          HyperedgeID& num_removed_single_pin_hyperedges,
                          HyperedgeVector& hyperedges,
                          parallel::scalable_vector<HyperedgeWeight>& hyperedges_weight,
                          parallel::scalable_vector<HypernodeWeight>& hypernodes_weight);

  Hypergraph readHypergraphFile(const std::string& filename,
                                const bool stable_construction_of_incident_edges = false);

  void readMetisFile(const std::string& filename,
                     HyperedgeID& num_hyperedges,
                     HypernodeID& num_hypernodes,
                     EdgeVector& hyperedges,
                     parallel::scalable_vector<HyperedgeWeight>& hyperedges_weight,
                     parallel::scalable_vector<HypernodeWeight>& hypernodes_weight);

  Hypergraph readMetisFile(const std::string& filename,
                           const bool stable_construction_of_incident_edges = false);

  Hypergraph readInputFile(const std::string& filename,
                           const FileFormat format,
                           const bool stable_construction_of_incident_edges = false);

  void readPartitionFile(const std::string& filename, std::vector<PartitionID>& partition);
  void writePartitionFile(const PartitionedHypergraph& phg, const std::string& filename);

}  // namespace io
}  // namespace mt_kahypar

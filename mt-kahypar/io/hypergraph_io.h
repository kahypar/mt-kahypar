/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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
                          parallel::scalable_vector<HypernodeWeight>& hypernodes_weight,
                          const bool remove_single_pin_hes = true);

  Hypergraph readHypergraphFile(const std::string& filename,
                                const bool stable_construction_of_incident_edges = false,
                                const bool remove_single_pin_hes = true);

  void readGraphFile(const std::string& filename,
                     HyperedgeID& num_hyperedges,
                     HypernodeID& num_hypernodes,
                     EdgeVector& hyperedges,
                     parallel::scalable_vector<HyperedgeWeight>& hyperedges_weight,
                     parallel::scalable_vector<HypernodeWeight>& hypernodes_weight);

  Hypergraph readGraphFile(const std::string& filename,
                           const bool stable_construction_of_incident_edges = false);

  Hypergraph readInputFile(const std::string& filename,
                           const FileFormat format,
                           const bool stable_construction_of_incident_edges = false,
                           const bool remove_single_pin_hes = true);

  void readPartitionFile(const std::string& filename, std::vector<PartitionID>& partition);
  void writePartitionFile(const PartitionedHypergraph& phg, const std::string& filename);

}  // namespace io
}  // namespace mt_kahypar

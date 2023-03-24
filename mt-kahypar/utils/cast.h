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

#include <string>

#include "include/libmtkahypartypes.h"

#include "mt-kahypar/macros.h"

namespace mt_kahypar::utils {

namespace {

std::string typeToString(const mt_kahypar_hypergraph_type_t type) {
  switch ( type ) {
    #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
    case STATIC_GRAPH: return "STATIC_GRAPH";
    #ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
    case DYNAMIC_GRAPH: return "DYNAMIC_GRAPH";
    #endif
    #endif
    case STATIC_HYPERGRAPH: return "STATIC_HYPERGRAPH";
    #ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
    case DYNAMIC_HYPERGRAPH: return "DYNAMIC_HYPERGRAPH";
    #endif
    case NULLPTR_HYPERGRAPH: return "NULLPTR_HYPERGRAPH";
  }
  return "UNDEFINED";
}

std::string typeToString(const mt_kahypar_partition_type_t type) {
  switch ( type ) {
    #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
    case MULTILEVEL_GRAPH_PARTITIONING: return "MULTILEVEL_GRAPH_PARTITIONING";
    #ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
    case N_LEVEL_GRAPH_PARTITIONING: return "N_LEVEL_GRAPH_PARTITIONING";
    #endif
    #endif
    case MULTILEVEL_HYPERGRAPH_PARTITIONING: return "MULTILEVEL_HYPERGRAPH_PARTITIONING";
    case LARGE_K_PARTITIONING: return "LARGE_K_PARTITIONING";
    #ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
    case N_LEVEL_HYPERGRAPH_PARTITIONING: return "N_LEVEL_HYPERGRAPH_PARTITIONING";
    #endif
    case NULLPTR_PARTITION: return "NULLPTR_PARTITION";
  }
  return "UNDEFINED";
}

} // namespace

template<typename Hypergraph>
Hypergraph& cast(mt_kahypar_hypergraph_t hypergraph) {
  if ( Hypergraph::TYPE != hypergraph.type ) {
    ERR("Cannot cast" << typeToString(hypergraph.type) << "to" << typeToString(Hypergraph::TYPE));
  }
  return *reinterpret_cast<Hypergraph*>(hypergraph.hypergraph);
}

template<typename Hypergraph>
const Hypergraph& cast_const(const mt_kahypar_hypergraph_t hypergraph) {
  if ( Hypergraph::TYPE != hypergraph.type ) {
    ERR("Cannot cast" << typeToString(hypergraph.type) << "to" << typeToString(Hypergraph::TYPE));
  }
  return *reinterpret_cast<const Hypergraph*>(hypergraph.hypergraph);
}

template<typename Hypergraph>
const Hypergraph& cast_const(mt_kahypar_hypergraph_const_t hypergraph) {
  if ( Hypergraph::TYPE != hypergraph.type ) {
    ERR("Cannot cast" << typeToString(hypergraph.type) << "to" << typeToString(Hypergraph::TYPE));
  }
  return *reinterpret_cast<const Hypergraph*>(hypergraph.hypergraph);
}

template<typename Hypergraph>
mt_kahypar_hypergraph_t hypergraph_cast(Hypergraph& hypergraph) {
  return mt_kahypar_hypergraph_t {
    reinterpret_cast<mt_kahypar_hypergraph_s*>(&hypergraph), Hypergraph::TYPE };
}

template<typename Hypergraph>
 mt_kahypar_hypergraph_const_t hypergraph_const_cast(const Hypergraph& hypergraph) {
  return mt_kahypar_hypergraph_const_t {
    reinterpret_cast<const mt_kahypar_hypergraph_s*>(&hypergraph), Hypergraph::TYPE };
}

template<typename PartitionedHypergraph>
PartitionedHypergraph& cast(mt_kahypar_partitioned_hypergraph_t phg) {
  if ( PartitionedHypergraph::TYPE != phg.type ) {
    ERR("Cannot cast" << typeToString(phg.type) << "to" << typeToString(PartitionedHypergraph::TYPE));
  }
  return *reinterpret_cast<PartitionedHypergraph*>(phg.partitioned_hg);
}

template<typename PartitionedHypergraph>
const PartitionedHypergraph& cast_const(const mt_kahypar_partitioned_hypergraph_t phg) {
  if ( PartitionedHypergraph::TYPE != phg.type ) {
    ERR("Cannot cast" << typeToString(phg.type) << "to" << typeToString(PartitionedHypergraph::TYPE));
  }
  return *reinterpret_cast<const PartitionedHypergraph*>(phg.partitioned_hg);
}

template<typename PartitionedHypergraph>
const PartitionedHypergraph& cast_const(mt_kahypar_partitioned_hypergraph_const_t phg) {
  if ( PartitionedHypergraph::TYPE != phg.type ) {
    ERR("Cannot cast" << typeToString(phg.type) << "to" << typeToString(PartitionedHypergraph::TYPE));
  }
  return *reinterpret_cast<const PartitionedHypergraph*>(phg.partitioned_hg);
}

template<typename PartitionedHypergraph>
mt_kahypar_partitioned_hypergraph_t partitioned_hg_cast(PartitionedHypergraph& phg) {
  return mt_kahypar_partitioned_hypergraph_t {
    reinterpret_cast<mt_kahypar_partitioned_hypergraph_s*>(&phg), PartitionedHypergraph::TYPE };
}

template<typename PartitionedHypergraph>
mt_kahypar_partitioned_hypergraph_const_t partitioned_hg_const_cast(const PartitionedHypergraph& phg) {
  return mt_kahypar_partitioned_hypergraph_const_t {
    reinterpret_cast<const mt_kahypar_partitioned_hypergraph_s*>(&phg), PartitionedHypergraph::TYPE };
}

}  // namespace mt_kahypar
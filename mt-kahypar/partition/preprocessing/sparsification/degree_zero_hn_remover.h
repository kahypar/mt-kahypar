/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <tbb/parallel_sort.h>

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/partition/refinement/rebalancing/rebalancer_common.h"
#include "mt-kahypar/weight/hypernode_weight_common.h"

namespace mt_kahypar {

template<typename TypeTraits>
class DegreeZeroHypernodeRemover {

  using Hypergraph = typename TypeTraits::Hypergraph;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;

 public:
  DegreeZeroHypernodeRemover(const Context& context) :
    _context(context),
    _removed_hns() { }

  DegreeZeroHypernodeRemover(const DegreeZeroHypernodeRemover&) = delete;
  DegreeZeroHypernodeRemover & operator= (const DegreeZeroHypernodeRemover &) = delete;

  DegreeZeroHypernodeRemover(DegreeZeroHypernodeRemover&&) = delete;
  DegreeZeroHypernodeRemover & operator= (DegreeZeroHypernodeRemover &&) = delete;

  // ! Remove all degree zero vertices
  HypernodeID removeDegreeZeroHypernodes(Hypergraph& hypergraph) {
    vec<double> weight_normalizer(hypergraph.dimension(), 0);
    for (Dimension d = 0; d < hypergraph.dimension(); ++d) {
      weight_normalizer[d] = 1 / static_cast<double>(hypergraph.totalWeight().at(d));
    }

    const HypernodeID current_num_nodes =
      hypergraph.initialNumNodes() - hypergraph.numRemovedHypernodes();
    HypernodeID num_removed_degree_zero_hypernodes = 0;
    for ( const HypernodeID& hn : hypergraph.nodes()  ) {
      if ( current_num_nodes - num_removed_degree_zero_hypernodes <= _context.coarsening.contraction_limit) {
        break;
      }
      if ( hypergraph.nodeDegree(hn) == 0 && !hypergraph.isFixed(hn) ) {
        hypergraph.removeDegreeZeroHypernode(hn);
        _removed_hns.emplace_back(hn, impl::normalizedSum(hypergraph.nodeWeight(hn), weight_normalizer));
        ++num_removed_degree_zero_hypernodes;
      }
    }
    return num_removed_degree_zero_hypernodes;
  }

  // ! Restore degree-zero vertices
  void restoreDegreeZeroHypernodes(PartitionedHypergraph& hypergraph) {
    // Sort degree-zero vertices in decreasing order of their weight
    tbb::parallel_sort(_removed_hns.begin(), _removed_hns.end(),
      [&](const std::pair<HypernodeID, float>& lhs, const std::pair<HypernodeID, float>& rhs) {
        return lhs.second > rhs.second || (lhs.second == rhs.second && lhs.first > rhs.first);
      });

    if (hypergraph.dimension() == 1) {
      // Sort blocks of partition in increasing order of their weight
      auto distance_to_max = [&](const PartitionID block) {
        return hypergraph.partWeight(block) - _context.partition.max_part_weights[block];
      };
      parallel::scalable_vector<PartitionID> blocks(_context.partition.k, 0);
      std::iota(blocks.begin(), blocks.end(), 0);
      std::sort(blocks.begin(), blocks.end(),
        [&](const PartitionID& lhs, const PartitionID& rhs) {
          return distance_to_max(lhs) < distance_to_max(rhs);
        });

      // Perform Bin-Packing
      for (const auto& [hn, _]: _removed_hns) {
        PartitionID to = blocks.front();
        hypergraph.restoreDegreeZeroHypernode(hn, to);
        PartitionID i = 0;
        while ( i + 1 < _context.partition.k &&
                distance_to_max(blocks[i]) > distance_to_max(blocks[i + 1]) ) {
          std::swap(blocks[i], blocks[i + 1]);
          ++i;
        }
      }
    } else {
      const vec<vec<double>> block_weight_normalizers = impl::computeBlockWeightNormalizers(_context);
      for (const auto& [hn, _]: _removed_hns) {
        const HNWeightConstRef weight = hypergraph.nodeWeight(hn);

        PartitionID best_target = kInvalidPartition;
        float best_block_rating = std::numeric_limits<float>::lowest();
        for (PartitionID to = 0; to < _context.partition.k; ++to) {
          float rating = impl::dotProduct(weight, _context.partition.max_part_weights[to] - hypergraph.partWeight(to), block_weight_normalizers[to]);
          if (!(hypergraph.partWeight(to) + weight <= _context.partition.max_part_weights[to])) {
            rating -= 1;
          }
          if (rating > best_block_rating) {
            best_target = to;
            best_block_rating = rating;
          }
        }
        ASSERT(best_target != kInvalidPartition);
        hypergraph.restoreDegreeZeroHypernode(hn, best_target);
      }
    }

    _removed_hns.clear();
  }

 private:
  const Context& _context;
  parallel::scalable_vector<std::pair<HypernodeID, float>> _removed_hns;
};

}  // namespace mt_kahypar

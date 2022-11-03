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

#include <vector>

#include "kahypar/utils/math.h"
#include "kahypar/utils/hash_vector.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/concurrent_bucket_map.h"
#include "mt-kahypar/datastructures/sparsifier_hypergraph.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/preprocessing/sparsification/i_hypergraph_sparsifier.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/utilities.h"

namespace mt_kahypar {
template<typename SimiliarNetCombiner>
class HypergraphSparsifier : public IHypergraphSparsifier {

  using Base = IHypergraphSparsifier;
  using SparsifierHypergraph = ds::SparsifierHypergraph<Hypergraph, HypergraphFactory>;

  using HashFunc = kahypar::math::MurmurHash<HypernodeID>;
  using HashValue = typename HashFunc::HashValue;
  using HashFuncVector = kahypar::HashFuncVector<HashFunc>;

  struct Footprint {
    explicit Footprint() :
      footprint(),
      he(kInvalidHyperedge) { }

    parallel::scalable_vector<HashValue> footprint;
    HyperedgeID he;

    bool operator==(const Footprint& other) const {
      ASSERT(footprint.size() == other.footprint.size());
      for ( size_t i = 0; i < footprint.size(); ++i ) {
        if ( footprint[i] != other.footprint[i] ) {
          return false;
        }
      }
      return true;
    }

    bool operator<(const Footprint& other) const {
      ASSERT(footprint.size() == other.footprint.size());
      for ( size_t i = 0; i < footprint.size(); ++i ) {
        if ( footprint[i] < other.footprint[i] ) {
          return true;
        } else if ( footprint[i] > other.footprint[i] ) {
          return false;
        }
      }
      return he < other.he;
    }

  };

  static constexpr bool enable_heavy_assert = false;

 public:
  HypergraphSparsifier(const Context& context) :
    Base(),
    _context(context),
    _sparsified_hg(),
    _sparsified_partitioned_hg(),
    _mapping() { }

  HypergraphSparsifier(const HypergraphSparsifier&) = delete;
  HypergraphSparsifier & operator= (const HypergraphSparsifier &) = delete;

  HypergraphSparsifier(HypergraphSparsifier&&) = delete;
  HypergraphSparsifier & operator= (HypergraphSparsifier &&) = delete;


 private:

  // ####################### Sparsification Functions #######################

  Hypergraph& sparsifiedHypergraphImpl() override final {
    return _sparsified_hg;
  }

  PartitionedHypergraph& sparsifiedPartitionedHypergraphImpl() override final {
    return _sparsified_partitioned_hg;
  }

  void sparsifyImpl(const Hypergraph& hypergraph) override final {
    ASSERT(_context.useSparsification());
    SparsifierHypergraph sparsified_hypergraph(hypergraph);
    utils::Timer& timer = utils::Utilities::instance().getTimer(_context.utility_id);

    // #################### STAGE 1 ####################
    // Heavy Hyperedge Removal
    // If the weight of all pins of a hyperedge is greater than a
    // certain threshold, we remove them from the hypergraph
    if ( _context.sparsification.use_heavy_net_removal ) {
      timer.start_timer("heavy_hyperedge_removal", "Heavy HE Removal");
      heavyHyperedgeRemovalSparsification(sparsified_hypergraph);
      timer.stop_timer("heavy_hyperedge_removal");
    }

    // #################### STAGE 2 ####################
    if ( _context.sparsification.use_similiar_net_removal ) {
      timer.start_timer("similiar_hyperedge_removal", "Similiar HE Removal");
      similiarHyperedgeRemoval(hypergraph, sparsified_hypergraph);
      timer.stop_timer("similiar_hyperedge_removal");
    }

    // #################### STAGE 3 ####################
    // Perform Degree-Zero Contractions
    // Degree-Zero hypernodes are contracted to supervertices such that
    // each supervertex has a weight smaller than the maximum allowed
    // node weight.
    if ( _context.sparsification.use_degree_zero_contractions ) {
      timer.start_timer("degree_zero_contraction", "Degree-Zero Contractions");
      degreeZeroSparsification(sparsified_hypergraph);
      timer.stop_timer("degree_zero_contraction");
    }

    // #################### STAGE 4 ####################
    // Construct sparsified hypergraph
    timer.start_timer("construct_sparsified_hypergraph", "Construct Sparsified HG");
    _sparsified_hg = sparsified_hypergraph.sparsify();
    _mapping = sparsified_hypergraph.getMapping();
    _sparsified_partitioned_hg = PartitionedHypergraph(
      _context.partition.k, _sparsified_hg, parallel_tag_t());
    timer.stop_timer("construct_sparsified_hypergraph");
  }

  void undoSparsificationImpl(PartitionedHypergraph& hypergraph) override final {
    hypergraph.doParallelForAllNodes([&](const HypernodeID hn) {
      ASSERT(hn < _mapping.size());
      const HypernodeID sparsified_hn = _mapping[hn];
      ASSERT(_sparsified_partitioned_hg.nodeIsEnabled(sparsified_hn));
      hypergraph.setNodePart(hn, _sparsified_partitioned_hg.partID(sparsified_hn));
    });
  }

 private:
  // ! Similiar to contractDegreeZeroHypernodes, but contraction is not applied
  // ! directly to hypergraph. Instead a mapping is computed that maps each vertex
  // ! of the original hypergraph to its supervertex and the weight of each
  // ! supervertex is aggregated in the hypernode weight vector.
  void degreeZeroSparsification(SparsifierHypergraph& hypergraph) {
    HypernodeID current_num_nodes = hypergraph.numNodes() - hypergraph.numRemovedNodes();
    HypernodeID degree_zero_supervertex = kInvalidHypernode;
    for (HypernodeID hn = 0; hn < hypergraph.numNodes(); ++hn) {
      if ( current_num_nodes <= _context.coarsening.contraction_limit ) {
        break;
      }

      if ( hypergraph.nodeDegree(hn) == 0 ) {
        bool was_removed = false;
        if ( degree_zero_supervertex != kInvalidHypernode ) {
          if ( hypergraph.nodeWeight(degree_zero_supervertex) +
               hypergraph.nodeWeight(hn) <=
               _context.coarsening.max_allowed_node_weight ) {
            // Remove vertex and aggregate its weight in its represenative supervertex
            hypergraph.contract(degree_zero_supervertex, hn);
            --current_num_nodes;
            was_removed = true;
          }
        }

        if ( !was_removed ) {
          degree_zero_supervertex = hn;
        }
      }
    }
  }

  // ! Removes hyperedges where the weight of all pins is greater
  // ! than a certain threshold. The threshold is specified in
  // ! '_context.initial_partitioning.max_hyperedge_pin_weight'
  void heavyHyperedgeRemovalSparsification(SparsifierHypergraph& hypergraph) {
    tbb::parallel_for(ID(0), hypergraph.numEdges(), [&](const HyperedgeID& e) {
      if ( hypergraph.edgeIsEnabled(e) ) {
        HypernodeWeight pin_weight = 0;
        for ( const HypernodeID& pin : hypergraph.pins(e) ) {
          pin_weight += hypergraph.nodeWeight(pin);
        }
        // Hyperedge will be include in sparsified hypergraph if its weight of
        // all pins is less than a predefined upper bound
        if ( pin_weight >= _context.sparsification.max_hyperedge_pin_weight ) {
          hypergraph.remove(e);
        }
      }
    });
  }

  void similiarHyperedgeRemoval(const Hypergraph& original_hg, SparsifierHypergraph& hypergraph) {
    HashFuncVector hash_functions(_context.sparsification.min_hash_footprint_size,
      utils::Randomize::instance().getRandomInt(0, 1000, sched_getcpu()));

    ds::ConcurrentBucketMap<Footprint> footprint_map;
    tbb::parallel_for(ID(0), hypergraph.numEdges(), [&](const HyperedgeID he) {
      if ( hypergraph.edgeIsEnabled(he) ) {
        Footprint he_footprint;
        he_footprint.footprint = {};
        he_footprint.he = he;
        for ( size_t i = 0; i < hash_functions.getHashNum(); ++i ) {
          he_footprint.footprint.push_back(minHash(hash_functions[i], hypergraph.pins(he)));
        }
        footprint_map.insert(combineHash(he_footprint), std::move(he_footprint));
      }
    });

    tbb::parallel_for(0UL, footprint_map.numBuckets(), [&](const size_t bucket) {
      auto& footprint_bucket = footprint_map.getBucket(bucket);
      if ( footprint_bucket.size() > 0 ) {
        std::sort(footprint_bucket.begin(), footprint_bucket.end());

        for ( size_t i = 0; i < footprint_bucket.size(); ++i ) {
          Footprint& representative = footprint_bucket[i];
          if ( representative.he != kInvalidHyperedge ) {
            parallel::scalable_vector<HypernodeID> rep_he = hypergraph.pins(representative.he);
            HyperedgeWeight rep_weight = hypergraph.edgeWeight(representative.he);
            bool exist_similiar_hes = false;
            for ( size_t j = i + 1; j < footprint_bucket.size(); ++j ) {
              Footprint& similiar_footprint = footprint_bucket[j];
              if ( similiar_footprint.he != kInvalidHyperedge ) {
                if ( representative == similiar_footprint ) {
                  const double jaccard_index = jaccard(
                    hypergraph.pins(representative.he), hypergraph.pins(similiar_footprint.he));
                  if ( jaccard_index >= _context.sparsification.jaccard_threshold ) {
                    rep_he = SimiliarNetCombiner::combine(original_hg, rep_he, hypergraph.pins(similiar_footprint.he));
                    rep_weight += hypergraph.edgeWeight(similiar_footprint.he);
                    hypergraph.remove(similiar_footprint.he);
                    similiar_footprint.he = kInvalidHyperedge;
                    exist_similiar_hes = true;
                  }
                } else {
                  break;
                }
              }
            }

            if ( exist_similiar_hes ) {
              hypergraph.replace(representative.he, std::move(rep_he));
              hypergraph.setEdgeWeight(representative.he, rep_weight);
            }
          }
        }
      }
      footprint_map.free(bucket);
    });
  }

  HashValue minHash(const HashFunc& hash_function,
                    const parallel::scalable_vector<HypernodeID>& hyperedge ) {
    HashValue hash_value = std::numeric_limits<HashValue>::max();
    for ( const HypernodeID& pin : hyperedge ) {
      hash_value = std::min(hash_value, hash_function(pin));
    }
    return hash_value;
  }

  HashValue combineHash(const Footprint& footprint) {
    HashValue hash_value = kEdgeHashSeed;
    for ( const HashValue& value : footprint.footprint ) {
      hash_value ^= value;
    }
    return hash_value;
  }

  double jaccard(const parallel::scalable_vector<HypernodeID>& lhs,
                 const parallel::scalable_vector<HypernodeID>& rhs) {
    const size_t min_size = std::min(lhs.size(), rhs.size());
    const size_t max_size = std::max(lhs.size(), rhs.size());
    if ( static_cast<double>(min_size) / static_cast<double>(max_size) <
         _context.sparsification.jaccard_threshold ) {
      return 0.0;
    }

    size_t intersection_size = 0;
    size_t i = 0;
    size_t j = 0;
    while ( i < lhs.size() && j < rhs.size() ) {
      if ( lhs[i] == rhs[j] ) {
        ++intersection_size;
        ++i;
        ++j;
      } else if ( lhs[i] < rhs[j] ) {
        ++i;
      } else {
        ++j;
      }
    }
    const size_t union_size = lhs.size() + rhs.size() - intersection_size;
    return static_cast<double>(intersection_size) /
      static_cast<double>(union_size);
  }

  const Context& _context;

  Hypergraph _sparsified_hg;
  PartitionedHypergraph _sparsified_partitioned_hg;
  parallel::scalable_vector<HypernodeID> _mapping;
};

}  // namespace mt_kahypar

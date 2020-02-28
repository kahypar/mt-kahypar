/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include <vector>

#include "kahypar/utils/math.h"
#include "kahypar/utils/hash_vector.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/clustering.h"
#include "mt-kahypar/datastructures/streaming_map.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
template<typename TypeTraits>
class HypergraphSparsifierT {

  using HyperGraph = typename TypeTraits::HyperGraph;
  using PartitionedHyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;
  using HyperGraphFactory = typename TypeTraits::HyperGraphFactory;
  using TBB = typename TypeTraits::TBB;
  using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;

  using HashFunc = kahypar::math::MurmurHash<HypernodeID>;
  using HashValue = typename HashFunc::HashValue;
  using HashFuncVector = kahypar::HashFuncVector<HashFunc>;

  struct Footprint {
    parallel::scalable_vector<HashValue> footprint;
    HyperedgeID he;

    bool operator==(const Footprint& other) {
      ASSERT(footprint.size() == other.footprint.size());
      for ( size_t i = 0; i < footprint.size(); ++i ) {
        if ( footprint[i] != other.footprint[i] ) {
          return false;
        }
      }
      return true;
    }

    bool operator<(const Footprint& other) {
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

  using FootprintMap = parallel::scalable_vector<parallel::scalable_vector<Footprint>>;

  static constexpr bool enable_heavy_assert = false;

 public:
  HypergraphSparsifierT(const Context& context,
                        const TaskGroupID task_group_id) :
    _context(context),
    _task_group_id(task_group_id),
    _is_init(false),
    _sparsified_hg(),
    _sparsified_partitioned_hg(),
    _removed_hes(),
    _removed_hns(),
    _mapping() { }

  HypergraphSparsifierT(const HypergraphSparsifierT&) = delete;
  HypergraphSparsifierT & operator= (const HypergraphSparsifierT &) = delete;

  HypergraphSparsifierT(HypergraphSparsifierT&&) = delete;
  HypergraphSparsifierT & operator= (HypergraphSparsifierT &&) = delete;


  // ####################### Sparsification Functions #######################

  bool isInitialized() const {
    return _is_init;
  }

  HyperGraph& sparsifiedHypergraph() {
    ASSERT(_is_init);
    return _sparsified_hg;
  }

  PartitionedHyperGraph& sparsifiedPartitionedHypergraph() {
    ASSERT(_is_init);
    return _sparsified_partitioned_hg;
  }

  void sparsify(const HyperGraph& hypergraph) {
    ASSERT(!_is_init);
    ASSERT(_mapping.size() == 0, "contractDegreeZeroHypernodes was triggered before");
    HyperedgeVector edge_vector;
    parallel::scalable_vector<int> vertices_to_numa_node;
    parallel::scalable_vector<HyperedgeWeight> hyperedge_weight;
    parallel::scalable_vector<HypernodeWeight> hypernode_weight;

    tbb::parallel_invoke([&] {
      _mapping.resize(hypergraph.initialNumNodes());
    }, [&]{
      hypernode_weight.resize(hypergraph.initialNumNodes());
    }, [&] {
      vertices_to_numa_node.resize(hypergraph.initialNumNodes());
    });

    // #################### STAGE 1 ####################
    // Perform Degree-Zero Contractions
    // Degree-Zero hypernodes are contracted to supervertices such that
    // each supervertex has a weight smaller than the maximum allowed
    // node weight.
    utils::Timer::instance().start_timer("degree_zero_contraction", "Degree-Zero Contractions");
    const HypernodeID num_hypernodes = degreeZeroSparsification(
      hypergraph, vertices_to_numa_node, hypernode_weight);
    utils::Timer::instance().stop_timer("degree_zero_contraction");

    // #################### STAGE 2 ####################
    // Heavy Hyperedge Removal
    // If the weight of all pins of a hyperedge is greater than a
    // certain threshold, we remove them from the hypergraph
    utils::Timer::instance().start_timer("heavy_hyperedge_removal", "Heavy HE Removal");
    HyperedgeID num_hyperedges = heavyHyperedgeRemovalSparsification(
      hypergraph, edge_vector, hyperedge_weight);
    utils::Timer::instance().stop_timer("heavy_hyperedge_removal");

    // #################### STAGE 3 ####################
    utils::Timer::instance().start_timer("similiar_hyperedge_removal", "Similiar HE Removal");
    similiarHyperedgeRemoval(edge_vector, hyperedge_weight);
    utils::Timer::instance().stop_timer("similiar_hyperedge_removal");

    // #################### STAGE 4 ####################
    // Construct sparsified hypergraph
    num_hyperedges = edge_vector.size();
    utils::Timer::instance().start_timer("construct_sparsified_hypergraph", "Construct Sparsified HG");
    if ( TBB::instance().num_used_numa_nodes() == 1 ) {
      _sparsified_hg = HyperGraphFactory::construct(
        _task_group_id, num_hypernodes, num_hyperedges,
        edge_vector, hyperedge_weight.data(), hypernode_weight.data());
    } else {
      _sparsified_hg = HyperGraphFactory::construct(
        _task_group_id, num_hypernodes, num_hyperedges,
        edge_vector, std::move(vertices_to_numa_node),
        hyperedge_weight.data(), hypernode_weight.data());
    }
    _sparsified_partitioned_hg = PartitionedHyperGraph(
      _context.partition.k, _task_group_id, _sparsified_hg);
    utils::Timer::instance().stop_timer("construct_sparsified_hypergraph");

    tbb::parallel_invoke([&] {
      parallel::parallel_free(edge_vector);
    }, [&] {
      parallel::parallel_free(hyperedge_weight, hypernode_weight);
    });

    _is_init = true;
  }


  // ! Removes single-pin hyperedges since they do not contribute
  // ! to cut or km1 metric.
  HyperedgeID removeSingleNodeHyperedges(HyperGraph& hypergraph) {
    ASSERT(!_is_init);
    HyperedgeID num_removed_single_node_hes = 0;
    for (const HyperedgeID& he : hypergraph.edges()) {
      if (hypergraph.edgeSize(he) == 1) {
        ++num_removed_single_node_hes;
        hypergraph.removeEdge(he);
        _removed_hes.push_back(he);
      }
    }
    return num_removed_single_node_hes;
  }

  // ! Contracts degree-zero vertices to degree-zero supervertices
  // ! We contract sets of degree-zero vertices such that the weight of
  // ! each supervertex is less than or equal than the maximum allowed
  // ! node weight for a vertex during coarsening.
  HypernodeID contractDegreeZeroHypernodes(HyperGraph& hypergraph) {
    ASSERT(!_is_init);
    _mapping.assign(hypergraph.initialNumNodes(), kInvalidHypernode);
    HypernodeID num_removed_degree_zero_hypernodes = 0;
    HypernodeID last_degree_zero_representative = kInvalidHypernode;
    HypernodeWeight last_degree_zero_weight = 0;
    for (const HypernodeID& hn : hypergraph.nodes()) {
      if ( hypergraph.nodeDegree(hn) == 0 ) {
        bool was_removed = false;
        if ( last_degree_zero_representative != kInvalidHypernode ) {
          const HypernodeWeight weight = hypergraph.nodeWeight(hn);
          if ( last_degree_zero_weight + weight <= _context.coarsening.max_allowed_node_weight ) {
            // Remove vertex and aggregate its weight in its represenative supervertex
            ++num_removed_degree_zero_hypernodes;
            hypergraph.removeHypernode(hn);
            _removed_hns.push_back(hn);
            was_removed = true;
            _mapping[hypergraph.originalNodeID(hn)] = last_degree_zero_representative;
            last_degree_zero_weight += weight;
            hypergraph.setNodeWeight(last_degree_zero_representative, last_degree_zero_weight);
          }
        }

        if ( !was_removed ) {
          last_degree_zero_representative = hn;
          last_degree_zero_weight = hypergraph.nodeWeight(hn);
        }
      }
    }
    return num_removed_degree_zero_hypernodes;
  }

  // ####################### Restore Operations #######################

  // ! Restore single-pin hyperedges
  void restoreSingleNodeHyperedges(PartitionedHyperGraph& hypergraph) {
    for (const HyperedgeID& he : _removed_hes) {
      hypergraph.restoreSinglePinHyperedge(he);
    }
  }

  // ! Restore degree-zero vertices
  // ! Each removed degree-zero vertex is assigned to the block of its supervertex.
  void restoreDegreeZeroHypernodes(PartitionedHyperGraph& hypergraph) {
    for ( const HypernodeID& hn : _removed_hns ) {
      const HypernodeID original_id = hypergraph.originalNodeID(hn);
      ASSERT(original_id < _mapping.size());
      const HypernodeID representative = _mapping[original_id];
      ASSERT(representative != kInvalidHypernode);
      // Restore degree-zero vertex and assign it to the block
      // of its supervertex
      hypergraph.enableHypernode(hn);
      hypergraph.setNodeWeight(representative,
        hypergraph.nodeWeight(representative) - hypergraph.nodeWeight(hn));
      hypergraph.setNodePart(hn, hypergraph.partID(representative));
    }
  }

  void undoSparsification(PartitionedHyperGraph& hypergraph) {
    ASSERT(_is_init);
    hypergraph.doParallelForAllNodes(_task_group_id, [&](const HypernodeID hn) {
      const HypernodeID original_id = hypergraph.originalNodeID(hn);
      ASSERT(original_id < _mapping.size());
      const HypernodeID original_sparsified_id = _mapping[original_id];
      const HypernodeID sparsified_hn = _sparsified_partitioned_hg.globalNodeID(original_sparsified_id);
      ASSERT(_sparsified_partitioned_hg.nodeIsEnabled(sparsified_hn));
      hypergraph.setNodePart(hn, _sparsified_partitioned_hg.partID(sparsified_hn));
    });
    hypergraph.initializeNumCutHyperedges(_task_group_id);
  }

  // ####################### Other #######################

  void assignAllDegreeZeroHypernodesToSameCommunity(HyperGraph& hypergraph, ds::Clustering& clustering) {
    ASSERT(hypergraph.initialNumNodes() == clustering.size());
    PartitionID community_id = -1;
    for ( const HypernodeID& hn : hypergraph.nodes() ) {
      if ( hypergraph.nodeDegree(hn) == 0 ) {
        if ( community_id >= 0 ) {
          clustering[hypergraph.originalNodeID(hn)] = community_id;
        } else {
          community_id = clustering[hypergraph.originalNodeID(hn)];
        }
      }
    }
  }

 private:
  // ! Similiar to contractDegreeZeroHypernodes, but contraction is not applied
  // ! directly to hypergraph. Instead a mapping is computed that maps each vertex
  // ! of the original hypergraph to its supervertex and the weight of each
  // ! supervertex is aggregated in the hypernode weight vector.
  HypernodeID degreeZeroSparsification(const HyperGraph& hypergraph,
                                       parallel::scalable_vector<int>& vertices_to_numa_node,
                                       parallel::scalable_vector<HypernodeWeight>& hypernode_weight) {
    ASSERT(hypergraph.initialNumNodes() == vertices_to_numa_node.size());
    ASSERT(hypergraph.initialNumNodes() == hypernode_weight.size());
    ASSERT(hypergraph.initialNumNodes() == _mapping.size());

    HypernodeID num_hypernodes = 0;
    HypernodeID last_degree_zero_representative = kInvalidHypernode;
    for (const HypernodeID& hn : hypergraph.nodes()) {
      const HypernodeID original_id = hypergraph.originalNodeID(hn);
      if ( hypergraph.nodeDegree(hn) == 0 ) {
        bool is_removed = false;
        if ( last_degree_zero_representative != kInvalidHypernode ) {
          const HypernodeWeight weight = hypergraph.nodeWeight(hn);
          if ( hypernode_weight[last_degree_zero_representative] + weight <=
               _context.coarsening.max_allowed_node_weight ) {
            // Aggregate weight of each degree-zero vertex in
            // its represenative supervertex
            _mapping[original_id] = last_degree_zero_representative;
            hypernode_weight[last_degree_zero_representative] += weight;
            is_removed = true;
          }
        }

        if ( !is_removed ) {
          hypernode_weight[num_hypernodes] = hypergraph.nodeWeight(hn);
          vertices_to_numa_node[num_hypernodes] = common::get_numa_node_of_vertex(hn);
          last_degree_zero_representative = num_hypernodes;
          _mapping[original_id] = num_hypernodes++;
        }
      } else {
        hypernode_weight[num_hypernodes] = hypergraph.nodeWeight(hn);
        vertices_to_numa_node[num_hypernodes] = common::get_numa_node_of_vertex(hn);
        _mapping[original_id] = num_hypernodes++;
      }
    }
    hypernode_weight.resize(num_hypernodes);
    vertices_to_numa_node.resize(num_hypernodes);

    return num_hypernodes;
  }

  // ! Removes hyperedges where the weight of all pins is greater
  // ! than a certain threshold. The threshold is specified in
  // ! '_context.initial_partitioning.max_hyperedge_pin_weight'
  HyperedgeID heavyHyperedgeRemovalSparsification(const HyperGraph& hypergraph,
                                                  HyperedgeVector& edge_vector,
                                                  parallel::scalable_vector<HyperedgeWeight>& hyperedge_weight) {
    ASSERT(hypergraph.initialNumNodes() == _mapping.size());

    utils::Timer::instance().start_timer("compute_included_hes", "Calc Included HEs");
    parallel::scalable_vector<HyperedgeID> include_in_sparsified_hg(hypergraph.initialNumEdges(), 0UL);
    hypergraph.doParallelForAllEdges(_task_group_id, [&](const HyperedgeID he) {
      HypernodeWeight pin_weight = 0;
      for ( const HypernodeID& pin : hypergraph.pins(he) ) {
        pin_weight += hypergraph.nodeWeight(pin);
      }
      // Hyperedge will be include in sparsified hypergraph if its weight of
      // all pins is less than a predefined upper bound
      if ( pin_weight <= _context.initial_partitioning.sparsification.max_hyperedge_pin_weight ) {
        include_in_sparsified_hg[hypergraph.originalEdgeID(he)] = 1UL;
      }
    });
    utils::Timer::instance().stop_timer("compute_included_hes");

    utils::Timer::instance().start_timer("sparsified_prefix_sum", "Sparsified HE Prefix Sum");
    parallel::TBBPrefixSum<HyperedgeID> sparsified_hg_prefix_sum(include_in_sparsified_hg);
    tbb::parallel_scan(tbb::blocked_range<size_t>(
      0UL, hypergraph.initialNumEdges()), sparsified_hg_prefix_sum);
    utils::Timer::instance().stop_timer("sparsified_prefix_sum");

    utils::Timer::instance().start_timer("add_sparsified_hes", "Add Sparsified HEs");
    tbb::parallel_invoke([&] {
      edge_vector.resize(sparsified_hg_prefix_sum.total_sum());
    }, [&]{
      hyperedge_weight.resize(sparsified_hg_prefix_sum.total_sum());
    });
    hypergraph.doParallelForAllEdges(_task_group_id, [&](const HyperedgeID he) {
      const HyperedgeID original_id = hypergraph.originalEdgeID(he);
      const bool included_in_sparsified_hg = sparsified_hg_prefix_sum.value(original_id);
      if ( included_in_sparsified_hg ) {
        const HyperedgeID pos = sparsified_hg_prefix_sum[original_id];
        hyperedge_weight[pos] = hypergraph.edgeWeight(he);
        for ( const HypernodeID& pin : hypergraph.pins(he) ) {
          edge_vector[pos].push_back(_mapping[hypergraph.originalNodeID(pin)]);
        }
        std::sort(edge_vector[pos].begin(), edge_vector[pos].end());
      }
    });
    utils::Timer::instance().stop_timer("add_sparsified_hes");

    return edge_vector.size();
  }

  void similiarHyperedgeRemoval(HyperedgeVector& edge_vector,
                                parallel::scalable_vector<HyperedgeWeight>& hyperedge_weight) {
    HashFuncVector hash_functions(_context.initial_partitioning.sparsification.min_hash_footprint_size,
      utils::Randomize::instance().getRandomInt(0, 1000, sched_getcpu()));

    ds::StreamingMap<HashValue, Footprint> hash_buckets;
    tbb::parallel_for(0UL, edge_vector.size(), [&](const HyperedgeID he) {
      Footprint he_footprint;
      he_footprint.footprint = {};
      he_footprint.he = he;
      for ( size_t i = 0; i < hash_functions.getHashNum(); ++i ) {
        he_footprint.footprint.push_back(minHash(hash_functions[i], edge_vector[he]));
      }
      hash_buckets.stream(combineHash(he_footprint), std::move(he_footprint));
    });

    FootprintMap footprint_map(hash_buckets.size());
    hash_buckets.copy(footprint_map, [&](const HashValue key) {
      return key % hash_buckets.size();
    });

    tbb::parallel_for(0UL, footprint_map.size(), [&](const size_t bucket) {
      parallel::scalable_vector<Footprint>& footprint_bucket = footprint_map[bucket];
      if ( footprint_bucket.size() > 0 ) {
        std::sort(footprint_bucket.begin(), footprint_bucket.end());

        for ( size_t i = 0; i < footprint_bucket.size(); ++i ) {
          Footprint& representative = footprint_bucket[i];
          if ( representative.he != kInvalidHyperedge ) {
            parallel::scalable_vector<HypernodeID> rep_he = edge_vector[representative.he];
            for ( size_t j = i + 1; j < footprint_bucket.size(); ++j ) {
              Footprint& similiar_footprint = footprint_bucket[j];
              if ( similiar_footprint.he != kInvalidHyperedge ) {
                if ( representative == similiar_footprint ) {
                  const double jaccard_index = jaccard(
                    edge_vector[representative.he], edge_vector[similiar_footprint.he]);
                  if ( jaccard_index >= _context.initial_partitioning.sparsification.jaccard_threshold ) {
                    sampleCombineHyperedges(rep_he, edge_vector[similiar_footprint.he]);
                    edge_vector[similiar_footprint.he].clear();
                    hyperedge_weight[representative.he] += hyperedge_weight[similiar_footprint.he];
                    hyperedge_weight[similiar_footprint.he] = 0;
                    similiar_footprint.he = kInvalidHyperedge;
                  }
                } else {
                  break;
                }
              }
            }
            edge_vector[representative.he] = std::move(rep_he);
          }
        }
      }
    });

    for ( size_t i = 0; i < edge_vector.size(); ++i ) {
      if ( edge_vector[i].size() == 0 ) {
        std::swap(edge_vector[i], edge_vector.back());
        std::swap(hyperedge_weight[i], hyperedge_weight.back());
        edge_vector.pop_back();
        hyperedge_weight.pop_back();
        --i;
      }
    }
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
         _context.initial_partitioning.sparsification.jaccard_threshold ) {
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

  void unionCombineHyperedges(parallel::scalable_vector<HypernodeID>& lhs,
                              const parallel::scalable_vector<HypernodeID>& rhs) {
    size_t lhs_size = lhs.size();
    size_t j = 0;
    for ( size_t i = 0; i < lhs_size; ++i ) {
      while ( j < rhs.size() && rhs[j] < lhs[i] ) {
        lhs.push_back(rhs[j++]);
      }
    }
    for ( ; j < rhs.size(); ++j ) {
      lhs.push_back(rhs[j]);
    }
    std::sort(lhs.begin(), lhs.end());
  }

  void sampleCombineHyperedges(parallel::scalable_vector<HypernodeID>& lhs,
                               const parallel::scalable_vector<HypernodeID>& rhs) {
    parallel::scalable_vector<HypernodeID> combined;
    int cpu_id = sched_getcpu();
    for ( const HypernodeID& pin : lhs ) {
      if ( utils::Randomize::instance().flipCoin(cpu_id) ) {
        combined.push_back(pin);
      }
    }
    for ( const HypernodeID& pin : rhs ) {
      if ( utils::Randomize::instance().flipCoin(cpu_id) ) {
        combined.push_back(pin);
      }
    }
    std::sort(combined.begin(), combined.end());
    combined.erase(std::unique(combined.begin(), combined.end()), combined.end());
    lhs = std::move(combined);
  }

  const Context& _context;
  const TaskGroupID _task_group_id;

  bool _is_init;
  HyperGraph _sparsified_hg;
  PartitionedHyperGraph _sparsified_partitioned_hg;
  parallel::scalable_vector<HyperedgeID> _removed_hes;
  parallel::scalable_vector<HypernodeID> _removed_hns;
  parallel::scalable_vector<HypernodeID> _mapping;
};
}  // namespace mt_kahypar

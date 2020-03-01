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

  class SparsifiedHypergraph {

    using AtomicNodeDegreeVector = parallel::scalable_vector<parallel::IntegralAtomicWrapper<HyperedgeID>>;
    using PinIterator = parallel::scalable_vector<HypernodeID>::const_iterator;

   public:
    SparsifiedHypergraph(const HyperGraph& hypergraph,
                          const TaskGroupID task_group_id) :
      _num_nodes(hypergraph.initialNumNodes()),
      _num_edges(hypergraph.initialNumEdges()),
      _task_group_id(task_group_id),
      _edge_vector(),
      _hyperedge_weight(),
      _he_included_in_sparsified_hg(),
      _mapping(),
      _vertices_to_numa_node(),
      _hypernode_weight(),
      _node_degrees(),
      _hn_included_in_sparsified_hg() {
      initialize(hypergraph, task_group_id);
    }

    SparsifiedHypergraph(const HypernodeID num_hypernodes,
                          const HyperedgeID num_hyperedges) :
      _num_nodes(num_hypernodes),
      _num_edges(num_hyperedges),
      _task_group_id(TBB::GLOBAL_TASK_GROUP),
      _edge_vector(),
      _hyperedge_weight(),
      _he_included_in_sparsified_hg(),
      _mapping(),
      _vertices_to_numa_node(),
      _hypernode_weight(),
      _node_degrees(),
      _hn_included_in_sparsified_hg() {
      tbb::parallel_invoke([&] {
        _edge_vector.resize(num_hyperedges);
      }, [&] {
        _hyperedge_weight.resize(num_hyperedges);
      }, [&] {
        _vertices_to_numa_node.resize(num_hypernodes);
      }, [&] {
        _hypernode_weight.resize(num_hypernodes);
      });
    }

    ~SparsifiedHypergraph() {
      tbb::parallel_invoke([&] {
        parallel::parallel_free(_edge_vector);
      }, [&] {
        parallel::parallel_free(_hyperedge_weight,
          _he_included_in_sparsified_hg, _vertices_to_numa_node,
          _hypernode_weight, _node_degrees,
          _hn_included_in_sparsified_hg);
      });
    }

    HypernodeID numNodes() const {
      return _num_nodes;
    }

    HyperedgeID numEdges() const {
      return _num_edges;
    }

    HypernodeWeight nodeWeight(const HypernodeID u) const {
      ASSERT(u < _num_nodes, V(u) << V(_num_nodes));
      return _hypernode_weight[u];
    }

    HyperedgeID nodeDegree(const HypernodeID u) const {
      ASSERT(u < _num_nodes);
      return _node_degrees[u];
    }

    bool nodeIsEnabled(const HypernodeID u) const {
      ASSERT(u < _num_nodes);
      return _hn_included_in_sparsified_hg[u];
    }

    void contract(const HypernodeID u, const HypernodeID v) {
      ASSERT(u != v);
      ASSERT(u < _num_nodes);
      ASSERT(v < _num_nodes);
      _mapping[v] = u;
      _hypernode_weight[u] += _hypernode_weight[v];
      _hypernode_weight[v] = 0;
      _hn_included_in_sparsified_hg[v] = 0;
    }

    parallel::scalable_vector<HypernodeID>&& getMapping() {
      return std::move(_mapping);
    }

    const parallel::scalable_vector<HypernodeID> pins(const HyperedgeID e) const {
      ASSERT(e < _num_edges);
      return _edge_vector[e];
    }

    HyperedgeWeight edgeWeight(const HyperedgeID e) const {
      ASSERT(e < _num_edges);
      return _hyperedge_weight[e];
    }

    void setEdgeWeight(const HyperedgeID e, const HyperedgeWeight weight) {
      ASSERT(e < _num_edges);
      _hyperedge_weight[e] = weight;
    }

    bool edgeIsEnabled(const HyperedgeID e) const {
      ASSERT(e < _num_edges);
      return _he_included_in_sparsified_hg[e];
    }

    void replace(const HyperedgeID e, parallel::scalable_vector<HyperedgeID>&& hyperedge) {
      ASSERT(e < _num_edges);
      for ( const HypernodeID& pin : _edge_vector[e]) {
        --_node_degrees[pin];
      }
      _edge_vector[e] = std::move(hyperedge);
      for ( const HypernodeID& pin : _edge_vector[e]) {
        ++_node_degrees[pin];
      }
    }

    void remove(const HyperedgeID e) {
      ASSERT(e < _num_edges);
      for ( const HypernodeID& pin : _edge_vector[e] ) {
        --_node_degrees[pin];
      }
      _edge_vector[e].clear();
      _hyperedge_weight[e] = 0;
      _he_included_in_sparsified_hg[e] = 0;
    }

    HyperGraph sparsify() {
      // Compute number of nodes and edges in sparsified hypergraph
      parallel::TBBPrefixSum<HyperedgeID> he_prefix_sum(_he_included_in_sparsified_hg);
      parallel::TBBPrefixSum<HyperedgeID> hn_prefix_sum(_hn_included_in_sparsified_hg);
      tbb::parallel_invoke([&] {
        tbb::parallel_scan(tbb::blocked_range<size_t>(ID(0), _num_edges), he_prefix_sum);
      }, [&] {
        tbb::parallel_scan(tbb::blocked_range<size_t>(ID(0), _num_nodes), hn_prefix_sum);
      });

      // Apply vertex id of sparsified hypergraph to mapping
      tbb::parallel_for(ID(0), _num_nodes, [&](const HypernodeID u) {
        _mapping[u] = hn_prefix_sum[_mapping[u]];
      });

      // Apply mapping to all enabled hyperedges
      tbb::parallel_for(ID(0), _num_edges, [&](const HyperedgeID& e) {
        if ( edgeIsEnabled(e) ) {
          parallel::scalable_vector<HypernodeID>& hyperedge = _edge_vector[e];
          for ( HypernodeID& pin : hyperedge ) {
            ASSERT(pin < _mapping.size());
            pin = _mapping[pin];
          }
          std::sort(hyperedge.begin(), hyperedge.end());
          hyperedge.erase(std::unique(hyperedge.begin(), hyperedge.end()), hyperedge.end());
        }
      });

      const HypernodeID num_hypernodes = hn_prefix_sum.total_sum();
      const HyperedgeID num_hyperedges = he_prefix_sum.total_sum();
      SparsifiedHypergraph sparsified_hypergraph(num_hypernodes, num_hyperedges);
      tbb::parallel_invoke([&] {
        tbb::parallel_for(ID(0), _num_edges, [&](const HyperedgeID e) {
          if ( he_prefix_sum.value(e) ) {
            const HyperedgeID sparsified_id = he_prefix_sum[e];
            sparsified_hypergraph._edge_vector[sparsified_id] = std::move(_edge_vector[e]);
            sparsified_hypergraph._hyperedge_weight[sparsified_id] = _hyperedge_weight[e];
          }
        });
      }, [&] {
        tbb::parallel_for(ID(0), _num_nodes, [&](const HypernodeID u) {
          if ( hn_prefix_sum.value(u) ) {
            const HyperedgeID sparsified_id = hn_prefix_sum[u];
            sparsified_hypergraph._vertices_to_numa_node[sparsified_id] = _vertices_to_numa_node[u];
            sparsified_hypergraph._hypernode_weight[sparsified_id] = _hypernode_weight[u];
          }
        });
      });

      if ( TBB::instance().num_used_numa_nodes() == 1 ) {
        return HyperGraphFactory::construct(
          _task_group_id, num_hypernodes, num_hyperedges,
          sparsified_hypergraph._edge_vector,
          sparsified_hypergraph._hyperedge_weight.data(),
          sparsified_hypergraph._hypernode_weight.data());
      } else {
        return HyperGraphFactory::construct(
          _task_group_id, num_hypernodes, num_hyperedges,
          sparsified_hypergraph._edge_vector,
          std::move(sparsified_hypergraph._vertices_to_numa_node),
          sparsified_hypergraph._hyperedge_weight.data(),
          sparsified_hypergraph._hypernode_weight.data());
      }
    }

   private:
    void initialize(const HyperGraph& hypergraph, const TaskGroupID task_group_id) {
      tbb::parallel_invoke([&] {
        tbb::parallel_invoke([&] {
          _edge_vector.resize(_num_edges);
        }, [&] {
          _hyperedge_weight.resize(_num_edges);
        }, [&] {
          _he_included_in_sparsified_hg.assign(_num_edges, 1);
        }, [&] {
          _node_degrees.assign(_num_nodes, parallel::IntegralAtomicWrapper<HyperedgeID>(0));
        });

        hypergraph.doParallelForAllEdges(task_group_id, [&](const HyperedgeID& he) {
          const HyperedgeID original_id = hypergraph.originalEdgeID(he);
          _hyperedge_weight[original_id] = hypergraph.edgeWeight(he);
          _edge_vector[original_id].reserve(hypergraph.edgeSize(he));
          for ( const HypernodeID& pin : hypergraph.pins(he) ) {
            const HypernodeID original_pin_id = hypergraph.originalNodeID(pin);
            ASSERT(original_pin_id < _num_nodes);
            _edge_vector[original_id].push_back(original_pin_id);
            ++_node_degrees[original_pin_id];
          }
          std::sort(_edge_vector[original_id].begin(), _edge_vector[original_id].end());
        });
      }, [&] {
        tbb::parallel_invoke([&] {
          _mapping.resize(_num_nodes);
        }, [&] {
          _vertices_to_numa_node.resize(_num_nodes);
        }, [&] {
          _hypernode_weight.resize(_num_nodes);
        }, [&] {
          _hn_included_in_sparsified_hg.assign(_num_nodes, 1);
        });

        hypergraph.doParallelForAllNodes(task_group_id, [&](const HypernodeID& hn) {
          const HypernodeID original_id = hypergraph.originalNodeID(hn);
          _mapping[original_id] = original_id;
          _vertices_to_numa_node[original_id] = common::get_numa_node_of_vertex(hn);
          _hypernode_weight[original_id] = hypergraph.nodeWeight(hn);
        });
      });
    }

    // Stats
    const HypernodeID _num_nodes;
    const HyperedgeID _num_edges;
    const TaskGroupID _task_group_id;

    // Hyperedges
    HyperedgeVector _edge_vector;
    parallel::scalable_vector<HyperedgeWeight> _hyperedge_weight;
    parallel::scalable_vector<HyperedgeID> _he_included_in_sparsified_hg;

    // Hypernodes
    parallel::scalable_vector<HypernodeID> _mapping;
    parallel::scalable_vector<int> _vertices_to_numa_node;
    parallel::scalable_vector<HypernodeWeight> _hypernode_weight;
    AtomicNodeDegreeVector _node_degrees;
    parallel::scalable_vector<HypernodeID> _hn_included_in_sparsified_hg;
  };


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
    SparsifiedHypergraph sparsified_hypergraph(hypergraph, _task_group_id);


    // #################### STAGE 1 ####################
    // Heavy Hyperedge Removal
    // If the weight of all pins of a hyperedge is greater than a
    // certain threshold, we remove them from the hypergraph
    utils::Timer::instance().start_timer("heavy_hyperedge_removal", "Heavy HE Removal");
    heavyHyperedgeRemovalSparsification(sparsified_hypergraph);
    utils::Timer::instance().stop_timer("heavy_hyperedge_removal");

    // #################### STAGE 2 ####################
    utils::Timer::instance().start_timer("similiar_hyperedge_removal", "Similiar HE Removal");
    similiarHyperedgeRemoval(sparsified_hypergraph);
    utils::Timer::instance().stop_timer("similiar_hyperedge_removal");

    // #################### STAGE 3 ####################
    // Perform Degree-Zero Contractions
    // Degree-Zero hypernodes are contracted to supervertices such that
    // each supervertex has a weight smaller than the maximum allowed
    // node weight.
    utils::Timer::instance().start_timer("degree_zero_contraction", "Degree-Zero Contractions");
    degreeZeroSparsification(sparsified_hypergraph);
    utils::Timer::instance().stop_timer("degree_zero_contraction");

    // #################### STAGE 4 ####################
    // Construct sparsified hypergraph
    utils::Timer::instance().start_timer("construct_sparsified_hypergraph", "Construct Sparsified HG");
    _sparsified_hg = sparsified_hypergraph.sparsify();
    _mapping = sparsified_hypergraph.getMapping();
    _sparsified_partitioned_hg = PartitionedHyperGraph(
      _context.partition.k, _task_group_id, _sparsified_hg);
    utils::Timer::instance().stop_timer("construct_sparsified_hypergraph");

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
  void degreeZeroSparsification(SparsifiedHypergraph& hypergraph) {
    ASSERT(!_is_init);
    HypernodeID degree_zero_supervertex = kInvalidHypernode;
    for (HypernodeID hn = 0; hn < hypergraph.numNodes(); ++hn) {
      if ( hypergraph.nodeDegree(hn) == 0 ) {
        bool was_removed = false;
        if ( degree_zero_supervertex != kInvalidHypernode ) {
          const HypernodeWeight weight = hypergraph.nodeWeight(hn);
          if ( hypergraph.nodeWeight(hn) + weight <= _context.coarsening.max_allowed_node_weight ) {
            // Remove vertex and aggregate its weight in its represenative supervertex
            hypergraph.contract(degree_zero_supervertex, hn);
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
  void heavyHyperedgeRemovalSparsification(SparsifiedHypergraph& hypergraph) {
    tbb::parallel_for(ID(0), hypergraph.numEdges(), [&](const HyperedgeID& e) {
      if ( hypergraph.edgeIsEnabled(e) ) {
        HypernodeWeight pin_weight = 0;
        for ( const HypernodeID& pin : hypergraph.pins(e) ) {
          pin_weight += hypergraph.nodeWeight(pin);
        }
        // Hyperedge will be include in sparsified hypergraph if its weight of
        // all pins is less than a predefined upper bound
        if ( pin_weight >= _context.initial_partitioning.sparsification.max_hyperedge_pin_weight ) {
          hypergraph.remove(e);
        }
      }
    });
  }

  void similiarHyperedgeRemoval(SparsifiedHypergraph& hypergraph) {
    HashFuncVector hash_functions(_context.initial_partitioning.sparsification.min_hash_footprint_size,
      utils::Randomize::instance().getRandomInt(0, 1000, sched_getcpu()));

    ds::StreamingMap<HashValue, Footprint> hash_buckets;
    tbb::parallel_for(ID(0), hypergraph.numEdges(), [&](const HyperedgeID he) {
      if ( hypergraph.edgeIsEnabled(he) ) {
        Footprint he_footprint;
        he_footprint.footprint = {};
        he_footprint.he = he;
        for ( size_t i = 0; i < hash_functions.getHashNum(); ++i ) {
          he_footprint.footprint.push_back(minHash(hash_functions[i], hypergraph.pins(he)));
        }
        hash_buckets.stream(combineHash(he_footprint), std::move(he_footprint));
      }
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
            parallel::scalable_vector<HypernodeID> rep_he = hypergraph.pins(representative.he);
            HyperedgeWeight rep_weight = hypergraph.edgeWeight(representative.he);
            bool exist_similiar_he = false;
            for ( size_t j = i + 1; j < footprint_bucket.size(); ++j ) {
              Footprint& similiar_footprint = footprint_bucket[j];
              if ( similiar_footprint.he != kInvalidHyperedge ) {
                if ( representative == similiar_footprint ) {
                  const double jaccard_index = jaccard(
                    hypergraph.pins(representative.he), hypergraph.pins(similiar_footprint.he));
                  if ( jaccard_index >= _context.initial_partitioning.sparsification.jaccard_threshold ) {
                    sampleCombineHyperedges(rep_he, hypergraph.pins(similiar_footprint.he));
                    rep_weight += hypergraph.edgeWeight(similiar_footprint.he);
                    hypergraph.remove(similiar_footprint.he);
                    similiar_footprint.he = kInvalidHyperedge;
                    exist_similiar_he = true;
                  }
                } else {
                  break;
                }
              }
            }

            if ( exist_similiar_he ) {
              hypergraph.replace(representative.he, std::move(rep_he));
              hypergraph.setEdgeWeight(representative.he, rep_weight);
            }
          }
        }
      }
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

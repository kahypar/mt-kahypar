/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Nikolai Maas <nikolai.maas@student.kit.edu>
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

#include "mt-kahypar/datastructures/dynamic_graph.h"

#include "tbb/blocked_range.h"
#include "tbb/parallel_scan.h"
#include "tbb/parallel_sort.h"
#include "tbb/parallel_reduce.h"
#include "tbb/concurrent_queue.h"

#include "mt-kahypar/parallel/stl/scalable_queue.h"
#include "mt-kahypar/datastructures/concurrent_bucket_map.h"
#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
namespace ds {

// ! Recomputes the total weight of the hypergraph (parallel)
void DynamicGraph::updateTotalWeight(parallel_tag_t) {
  _total_weight = tbb::parallel_reduce(tbb::blocked_range<HypernodeID>(ID(0), _num_nodes), 0,
    [this](const tbb::blocked_range<HypernodeID>& range, HypernodeWeight init) {
      HypernodeWeight weight = init;
      for (HypernodeID hn = range.begin(); hn < range.end(); ++hn) {
        if ( nodeIsEnabled(hn) ) {
          weight += this->_nodes[hn].weight();
        }
      }
      return weight;
    }, std::plus<HypernodeWeight>()) + _removed_degree_zero_hn_weight;
}

// ! Recomputes the total weight of the hypergraph (sequential)
void DynamicGraph::updateTotalWeight() {
  _total_weight = 0;
  for ( const HypernodeID& hn : nodes() ) {
    if ( nodeIsEnabled(hn) ) {
      _total_weight += nodeWeight(hn);
    }
  }
  _total_weight += _removed_degree_zero_hn_weight;
}

/**!
 * Registers a contraction in the hypergraph whereas vertex u is the representative
 * of the contraction and v its contraction partner. Several threads can call this function
 * in parallel. The function adds the contraction of u and v to a contraction tree that determines
 * a parallel execution order and synchronization points for all running contractions.
 * The contraction can be executed by calling function contract(v, max_node_weight).
 */
bool DynamicGraph::registerContraction(const HypernodeID u, const HypernodeID v) {
  return _contraction_tree.registerContraction(u, v, _version,
                                               [&](HypernodeID u) { acquireHypernode(u); },
                                               [&](HypernodeID u) { releaseHypernode(u); });
}

/**!
 * Contracts a previously registered contraction. Representative u of vertex v is looked up
 * in the contraction tree and performed if there are no pending contractions in the subtree
 * of v and the contractions respects the maximum allowed node weight. If (u,v) is the last
 * pending contraction in the subtree of u then the function recursively contracts also
 * u (if any contraction is registered). Therefore, function can return several contractions
 * or also return an empty contraction vector.
 */
size_t DynamicGraph::contract(const HypernodeID v,
                              const HypernodeWeight max_node_weight) {
  ASSERT(_contraction_tree.parent(v) != v, "No contraction registered for node " << v);

  HypernodeID x = _contraction_tree.parent(v);
  HypernodeID y = v;
  ContractionResult res = ContractionResult::CONTRACTED;
  size_t num_contractions = 0;
  // We perform all contractions registered in the contraction tree
  // as long as there are no pending contractions
  while ( x != y && res != ContractionResult::PENDING_CONTRACTIONS) {
    // Perform Contraction
    res = contract(x, y, max_node_weight);
    if ( res == ContractionResult::CONTRACTED ) {
      ++num_contractions;
    }
    y = x;
    x = _contraction_tree.parent(y);
  }
  return num_contractions;
}


 /**
   * Uncontracts a batch of contractions in parallel. The batches must be uncontracted exactly
   * in the order computed by the function createBatchUncontractionHierarchy(...).
   * The two uncontraction functions are required by the partitioned graph to update
   * gain cache values.
   */
void DynamicGraph::uncontract(const Batch& batch,
                                   const UncontractionFunction& case_one_func,
                                   const UncontractionFunction& case_two_func) {
  ASSERT(batch.size() > 0UL);
  ASSERT([&] {
    const HypernodeID expected_batch_index = hypernode(batch[0].v).batchIndex();
    for ( const Memento& memento : batch ) {
      if ( hypernode(memento.v).batchIndex() != expected_batch_index ) {
        LOG << "Batch contains uncontraction from different batches."
            << "Hypernode" << memento.v << "with version" << hypernode(memento.v).batchIndex()
            << "but expected is" << expected_batch_index;
        return false;
      }
      if ( _contraction_tree.version(memento.v) != _version ) {
        LOG << "Batch contains uncontraction from a different version."
            << "Hypernode" << memento.v << "with version" << _contraction_tree.version(memento.v)
            << "but expected is" << _version;
        return false;
      }
    }
    return true;
  }(), "Batch contains uncontractions from different batches or from a different hypergraph version");

  tbb::parallel_for(0UL, batch.size(), [&](const size_t i) {
    const Memento& memento = batch[i];
    ASSERT(!hypernode(memento.u).isDisabled(), "Hypernode" << memento.u << "is disabled");
    ASSERT(hypernode(memento.v).isDisabled(), "Hypernode" << memento.v << "is not invalid");

    // Restore incident net list of u and v
    _adjacency_array.uncontract(memento.u, memento.v,
      [&](const HyperedgeID e) {
        case_one_func(memento.u, memento.v, e);
      }, [&](const HyperedgeID e) {
        case_two_func(memento.u, memento.v, e);
      }, [&](const HypernodeID u) {
        acquireHypernode(u);
      }, [&](const HypernodeID u) {
        releaseHypernode(u);
      });

    acquireHypernode(memento.u);
    // Restore hypernode v which includes enabling it and subtract its weight
    // from its representative
    hypernode(memento.v).enable();
    hypernode(memento.u).setWeight(hypernode(memento.u).weight() - hypernode(memento.v).weight());
    releaseHypernode(memento.u);
  });
}

/**
 * Computes a batch uncontraction hierarchy. A batch is a vector of mementos
 * (uncontractions) that are uncontracted in parallel. The function returns a vector
 * of versioned batch vectors. A new version of the hypergraph is induced if we perform
 * single-pin and parallel net detection. Once we process all batches of a versioned
 * batch vector, we have to restore all previously removed single-pin and parallel nets
 * in order to uncontract the next batch vector. We create for each version of the
 * hypergraph a seperate batch uncontraction hierarchy (see createBatchUncontractionHierarchyOfVersion(...))
 */
VersionedBatchVector DynamicGraph::createBatchUncontractionHierarchy(const size_t batch_size) {
  const size_t num_versions = _version + 1;
  utils::Timer::instance().start_timer("finalize_contraction_tree", "Finalize Contraction Tree");
  // Finalizes the contraction tree such that it is traversable in a top-down fashion
  // and contains subtree size for each  tree node
  _contraction_tree.finalize(num_versions);
  utils::Timer::instance().stop_timer("finalize_contraction_tree");

  utils::Timer::instance().start_timer("create_versioned_batches", "Create Versioned Batches");
  VersionedBatchVector versioned_batches(num_versions);
  parallel::scalable_vector<size_t> batch_sizes_prefix_sum(num_versions, 0);
  BatchIndexAssigner batch_index_assigner(_num_nodes, batch_size);
  for ( size_t version = 0; version < num_versions; ++version ) {
    versioned_batches[version] =
      _contraction_tree.createBatchUncontractionHierarchyForVersion(batch_index_assigner, version);
    if ( version > 0 ) {
      batch_sizes_prefix_sum[version] =
        batch_sizes_prefix_sum[version - 1] + versioned_batches[version - 1].size();
    }
    batch_index_assigner.reset(versioned_batches[version].size());
  }
  utils::Timer::instance().stop_timer("create_versioned_batches");

  return versioned_batches;
}

/**
 * Removes single-pin and parallel nets from the hypergraph. The weight
 * of a set of identical nets is aggregated in one representative hyperedge
 * and single-pin hyperedges are removed. Returns a vector of removed hyperedges.
 */
parallel::scalable_vector<DynamicGraph::ParallelHyperedge> DynamicGraph::removeSinglePinAndParallelHyperedges() {
  ++_version;
  return _adjacency_array.removeSinglePinAndParallelEdges();
}

/**
 * Restores a previously removed set of singple-pin and parallel hyperedges. Note, that hes_to_restore
 * must be exactly the same and given in the reverse order as returned by removeSinglePinAndParallelNets(...).
 */
void DynamicGraph::restoreSinglePinAndParallelNets(const parallel::scalable_vector<ParallelHyperedge>& hes_to_restore) {
  _adjacency_array.restoreSinglePinAndParallelEdges(hes_to_restore);
  --_version;
}

// ! Copy dynamic hypergraph in parallel
DynamicGraph DynamicGraph::copy(parallel_tag_t) {
  DynamicGraph hypergraph;

  hypergraph._num_nodes = _num_nodes;
  hypergraph._num_removed_nodes = _num_removed_nodes;
  hypergraph._removed_degree_zero_hn_weight = _removed_degree_zero_hn_weight;
  hypergraph._num_edges = _num_edges;
  hypergraph._total_weight = _total_weight;
  hypergraph._version = _version;
  hypergraph._contraction_index.store(_contraction_index.load());

  tbb::parallel_invoke([&] {
    hypergraph._nodes.resize(_nodes.size());
    memcpy(hypergraph._nodes.data(), _nodes.data(),
      sizeof(Hypernode) * _nodes.size());
  }, [&] {
    hypergraph._adjacency_array = _adjacency_array.copy(parallel_tag_t());
  }, [&] {
    hypergraph._acquired_nodes.resize(_acquired_nodes.size());
    tbb::parallel_for(ID(0), _num_nodes, [&](const HypernodeID& hn) {
      hypergraph._acquired_nodes[hn] = _acquired_nodes[hn];
    });
  }, [&] {
    hypergraph._contraction_tree = _contraction_tree.copy(parallel_tag_t());
  });
  return hypergraph;
}

// ! Copy dynamic hypergraph sequential
DynamicGraph DynamicGraph::copy() {
  DynamicGraph hypergraph;

  hypergraph._num_nodes = _num_nodes;
  hypergraph._num_removed_nodes = _num_removed_nodes;
  hypergraph._removed_degree_zero_hn_weight = _removed_degree_zero_hn_weight;
  hypergraph._num_edges = _num_edges;
  hypergraph._total_weight = _total_weight;
  hypergraph._version = _version;
  hypergraph._contraction_index.store(_contraction_index.load());

  hypergraph._nodes.resize(_nodes.size());
  memcpy(hypergraph._nodes.data(), _nodes.data(),
    sizeof(Hypernode) * _nodes.size());
    hypergraph._adjacency_array = _adjacency_array.copy(parallel_tag_t());
  hypergraph._acquired_nodes.resize(_num_nodes);
  for ( HypernodeID hn = 0; hn < _num_nodes; ++hn ) {
    hypergraph._acquired_nodes[hn] = _acquired_nodes[hn];
  }
  hypergraph._contraction_tree = _contraction_tree.copy();

  return hypergraph;
}

void DynamicGraph::memoryConsumption(utils::MemoryTreeNode* parent) const {
  ASSERT(parent);

  parent->addChild("Hypernodes", sizeof(Hypernode) * _nodes.size());
  parent->addChild("Incident Nets", _adjacency_array.size_in_bytes());
  parent->addChild("Hypernode Ownership Vector", sizeof(bool) * _acquired_nodes.size());

  utils::MemoryTreeNode* contraction_tree_node = parent->addChild("Contraction Tree");
  _contraction_tree.memoryConsumption(contraction_tree_node);
}

// ! Only for testing
// bool DynamicGraph::verifyIncidenceArrayAndIncidentNets() {
//   bool success = true;
//   tbb::parallel_invoke([&] {
//     doParallelForAllNodes([&](const HypernodeID& hn) {
//       for ( const HyperedgeID& he : incidentEdges(hn) ) {
//         bool found = false;
//         for ( const HypernodeID& pin : pins(he) ) {
//           if ( pin == hn ) {
//             found = true;
//             break;
//           }
//         }
//         if ( !found ) {
//           LOG << "Hypernode" << hn << "not found in incidence array of net" << he;
//           success = false;
//         }
//       }
//     });
//   }, [&] {
//     doParallelForAllEdges([&](const HyperedgeID& he) {
//       for ( const HypernodeID& pin : pins(he) ) {
//         bool found = false;
//         for ( const HyperedgeID& e : incidentEdges(pin) ) {
//           if ( e == he ) {
//             found = true;
//             break;
//           }
//         }
//         if ( !found ) {
//           LOG << "Hyperedge" << he << "not found in incident nets of vertex" << pin;
//           success = false;
//         }
//       }
//     });
//   });
//   return success;
// }

/**!
 * Contracts a previously registered contraction. The contraction of u and v is
 * performed if there are no pending contractions in the subtree of v and the
 * contractions respects the maximum allowed node weight. In case the contraction
 * was performed successfully, enum type CONTRACTED is returned. If contraction
 * was not performed either WEIGHT_LIMIT_REACHED (in case sum of both vertices is
 * greater than the maximum allowed node weight) or PENDING_CONTRACTIONS (in case
 * there are some unfinished contractions in the subtree of v) is returned.
 */
DynamicGraph::ContractionResult DynamicGraph::contract(const HypernodeID u,
                                                       const HypernodeID v,
                                                       const HypernodeWeight max_node_weight) {

  // Acquire ownership in correct order to prevent deadlocks
  if ( u < v ) {
    acquireHypernode(u);
    acquireHypernode(v);
  } else {
    acquireHypernode(v);
    acquireHypernode(u);
  }

  // Contraction is valid if
  //  1.) Contraction partner v is enabled
  //  2.) There are no pending contractions on v
  //  4.) Resulting node weight is less or equal than a predefined upper bound
  const bool contraction_partner_valid =
    nodeIsEnabled(v) && _contraction_tree.pendingContractions(v) == 0;
  const bool less_or_equal_than_max_node_weight =
    hypernode(u).weight() + hypernode(v).weight() <= max_node_weight;
  if ( contraction_partner_valid && less_or_equal_than_max_node_weight ) {
    ASSERT(nodeIsEnabled(u), "Hypernode" << u << "is disabled!");
    hypernode(u).setWeight(nodeWeight(u) + nodeWeight(v));
    hypernode(v).disable();
    releaseHypernode(u);
    releaseHypernode(v);

    HypernodeID contraction_start = _contraction_index.load();

    // Contract incident net lists of u and v
    _adjacency_array.contract(u, v, [&](const HypernodeID u) {
      acquireHypernode(u);
    }, [&](const HypernodeID u) {
      releaseHypernode(u);
    });

    HypernodeID contraction_end = ++_contraction_index;
    acquireHypernode(u);
    _contraction_tree.unregisterContraction(u, v, contraction_start, contraction_end);
    releaseHypernode(u);
    return ContractionResult::CONTRACTED;
  } else {
    ContractionResult res = ContractionResult::PENDING_CONTRACTIONS;
    if ( !less_or_equal_than_max_node_weight && nodeIsEnabled(v) &&
         _contraction_tree.parent(v) == u ) {
      _contraction_tree.unregisterContraction(u, v,
        kInvalidHypernode, kInvalidHypernode, true /* failed */);
      res = ContractionResult::WEIGHT_LIMIT_REACHED;
    }
    releaseHypernode(u);
    releaseHypernode(v);
    return res;
  }
}

// bool DynamicGraph::verifyBatchIndexAssignments(
//   const BatchIndexAssigner& batch_assigner,
//   const parallel::scalable_vector<parallel::scalable_vector<BatchAssignment>>& local_batch_assignments) const {
//   parallel::scalable_vector<BatchAssignment> assignments;
//   for ( size_t i = 0; i < local_batch_assignments.size(); ++i ) {
//     for ( const BatchAssignment& batch_assign : local_batch_assignments[i] ) {
//       assignments.push_back(batch_assign);
//     }
//   }
//   std::sort(assignments.begin(), assignments.end(),
//     [&](const BatchAssignment& lhs, const BatchAssignment& rhs) {
//       return lhs.batch_index < rhs.batch_index ||
//         (lhs.batch_index == rhs.batch_index && lhs.batch_pos < rhs.batch_pos);
//     });

//   if ( assignments.size() > 0 ) {
//     if ( assignments[0].batch_index != 0 || assignments[0].batch_pos != 0 ) {
//       LOG << "First uncontraction should start at batch 0 at position 0"
//           << V(assignments[0].batch_index) << V(assignments[0].batch_pos);
//       return false;
//     }

//     for ( size_t i = 1; i < assignments.size(); ++i ) {
//       if ( assignments[i - 1].batch_index == assignments[i].batch_index ) {
//         if ( assignments[i - 1].batch_pos + 1 != assignments[i].batch_pos ) {
//           LOG << "Batch positions are not consecutive"
//               << V(i) << V(assignments[i - 1].batch_pos) << V(assignments[i].batch_pos);
//           return false;
//         }
//       } else {
//         if ( assignments[i - 1].batch_index + 1 != assignments[i].batch_index ) {
//           LOG << "Batch indices are not consecutive"
//               << V(i) << V(assignments[i - 1].batch_index) << V(assignments[i].batch_index);
//           return false;
//         }
//         if ( assignments[i].batch_pos != 0 ) {
//           LOG << "First uncontraction of each batch should start at position 0"
//               << V(assignments[i].batch_pos);
//           return false;
//         }
//         if ( assignments[i - 1].batch_pos + 1 != batch_assigner.batchSize(assignments[i - 1].batch_index) ) {
//           LOG << "Position of last uncontraction in batch" << assignments[i - 1].batch_index
//               << "does not match size of batch"
//               << V(assignments[i - 1].batch_pos) << V(batch_assigner.batchSize(assignments[i - 1].batch_index));
//           return false;
//         }
//       }
//     }
//   }

//   return true;
// }

} // namespace ds
} // namespace mt_kahypar
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
  // Acquires ownership of vertex v that gives the calling thread exclusive rights
  // to modify the contraction tree entry of v
  acquireHypernode(v);

  // If there is no other contraction registered for vertex v
  // we try to determine its representative in the contraction tree
  if ( _contraction_tree.parent(v) == v ) {

    HypernodeID w = u;
    bool cycle_detected = false;
    while ( true ) {
      // Search for representative of u in the contraction tree.
      // It is either a root of the contraction tree or a vertex
      // with a reference count greater than zero, which indicates
      // that there are still ongoing contractions on this node that
      // have to be processed.
      while ( _contraction_tree.parent(w) != w &&
              _contraction_tree.pendingContractions(w) == 0 ) {
        w = _contraction_tree.parent(w);
        if ( w == v ) {
          cycle_detected = true;
          break;
        }
      }

      if ( !cycle_detected ) {
        // In case contraction of u and v does not induce any
        // cycle in the contraction tree we try to acquire vertex w
        if ( w < v ) {
          // Acquire ownership in correct order to prevent deadlocks
          releaseHypernode(v);
          acquireHypernode(w);
          acquireHypernode(v);
          if ( _contraction_tree.parent(v) != v ) {
            releaseHypernode(v);
            releaseHypernode(w);
            return false;
          }
        } else {
          acquireHypernode(w);
        }

        // Double-check condition of while loop above after acquiring
        // ownership of w
        if ( _contraction_tree.parent(w) != w &&
              _contraction_tree.pendingContractions(w) == 0 ) {
          // In case something changed, we release ownership of w and
          // search again for the representative of u.
          releaseHypernode(w);
        } else {
          // Otherwise we perform final cycle check to verify that
          // contraction of u and v will not introduce any new cycle.
          HypernodeID x = w;
          do {
            x = _contraction_tree.parent(x);
            if ( x == v ) {
              cycle_detected = true;
              break;
            }
          } while ( _contraction_tree.parent(x) != x );

          if ( cycle_detected ) {
            releaseHypernode(w);
            releaseHypernode(v);
            return false;
          }

          // All checks succeded, we can safely increment the
          // reference count of w and update the contraction tree
          break;
        }
      } else {
        releaseHypernode(v);
        return false;
      }
    }

    // Increment reference count of w indicating that there pending
    // contraction at vertex w and update contraction tree.
    _contraction_tree.registerContraction(w, v, _version);

    releaseHypernode(w);
    releaseHypernode(v);
    return true;
  } else {
    releaseHypernode(v);
    return false;
  }
}

/**!
 * Contracts a previously registered contraction. Representative u of vertex v is looked up
 * in the contraction tree and performed if there are no pending contractions in the subtree
 * of v and the contractions respects the maximum allowed node weight. If (u,v) is the last
 * pending contraction in the subtree of u then the function recursively contracts also
 * u (if any contraction is registered). Therefore, function can return several contractions
 * or also return an empty contraction vector.
 */
// size_t DynamicGraph::contract(const HypernodeID v,
//                                    const HypernodeWeight max_node_weight) {
//   ASSERT(_contraction_tree.parent(v) != v, "No contraction registered for hypernode" << v);

//   HypernodeID x = _contraction_tree.parent(v);
//   HypernodeID y = v;
//   ContractionResult res = ContractionResult::CONTRACTED;
//   size_t num_contractions = 0;
//   // We perform all contractions registered in the contraction tree
//   // as long as there are no pending contractions
//   while ( x != y && res != ContractionResult::PENDING_CONTRACTIONS) {
//     // Perform Contraction
//     res = contract(x, y, max_node_weight);
//     if ( res == ContractionResult::CONTRACTED ) {
//       ++num_contractions;
//     }
//     y = x;
//     x = _contraction_tree.parent(y);
//   }
//   return num_contractions;
// }


/**
 * Uncontracts a batch of contractions in parallel. The batches must be uncontracted exactly
 * in the order computed by the function createBatchUncontractionHierarchy(...).
 * The two uncontraction functions are required by the partitioned hypergraph to restore
 * pin counts and gain cache values.
 */
// void DynamicGraph::uncontract(const Batch& batch,
//                                    const UncontractionFunction& case_one_func,
//                                    const UncontractionFunction& case_two_func) {
//   ASSERT(batch.size() > 0UL);
//   ASSERT([&] {
//     const HypernodeID expected_batch_index = hypernode(batch[0].v).batchIndex();
//     for ( const Memento& memento : batch ) {
//       if ( hypernode(memento.v).batchIndex() != expected_batch_index ) {
//         LOG << "Batch contains uncontraction from different batches."
//             << "Hypernode" << memento.v << "with version" << hypernode(memento.v).batchIndex()
//             << "but expected is" << expected_batch_index;
//         return false;
//       }
//       if ( _contraction_tree.version(memento.v) != _version ) {
//         LOG << "Batch contains uncontraction from a different version."
//             << "Hypernode" << memento.v << "with version" << _contraction_tree.version(memento.v)
//             << "but expected is" << _version;
//         return false;
//       }
//     }
//     return true;
//   }(), "Batch contains uncontractions from different batches or from a different hypergraph version");

//   _hes_to_resize_flag_array.reset();
//   tbb::parallel_for(0UL, batch.size(), [&](const size_t i) {
//     const Memento& memento = batch[i];
//     ASSERT(!hypernode(memento.u).isDisabled(), "Hypernode" << memento.u << "is disabled");
//     ASSERT(hypernode(memento.v).isDisabled(), "Hypernode" << memento.v << "is not invalid");

//     // Restore incident net list of u and v
//     const HypernodeID batch_index = hypernode(batch[0].v).batchIndex();
//     _incident_nets.uncontract(memento.u, memento.v,
//       [&](const HyperedgeID e) {
//         // In that case, u and v were both previously part of hyperedge e.
//         if ( !_hes_to_resize_flag_array[e] &&
//              _hes_to_resize_flag_array.compare_and_set_to_true(e) ) {
//           // This part is only triggered once for each hyperedge per batch uncontraction.
//           // It restores all pins that are part of the current batch as contraction partners
//           // in hyperedge e
//           restoreHyperedgeSizeForBatch(e, batch_index);
//         }

//         acquireEdge(e);
//         case_one_func(memento.u, memento.v, e);
//         releaseEdge(e);
//       }, [&](const HyperedgeID e) {
//         // In that case only v was part of hyperedge e before and
//         // u must be replaced by v in hyperedge e
//         const size_t slot_of_u = findPositionOfPinInIncidenceArray(memento.u, e);

//         acquireEdge(e);
//         ASSERT(_incidence_array[slot_of_u] == memento.u);
//         _incidence_array[slot_of_u] = memento.v;
//         case_two_func(memento.u, memento.v, e);
//         releaseEdge(e);
//       }, [&](const HypernodeID u) {
//         acquireHypernode(u);
//       }, [&](const HypernodeID u) {
//         releaseHypernode(u);
//       });

//     acquireHypernode(memento.u);
//     // Restore hypernode v which includes enabling it and subtract its weight
//     // from its representative
//     hypernode(memento.v).enable();
//     hypernode(memento.u).setWeight(hypernode(memento.u).weight() - hypernode(memento.v).weight());
//     releaseHypernode(memento.u);
//   });
// }

/**
 * Computes a batch uncontraction hierarchy. A batch is a vector of mementos
 * (uncontractions) that are uncontracted in parallel. The function returns a vector
 * of versioned batch vectors. A new version of the hypergraph is induced if we perform
 * single-pin and parallel net detection. Once we process all batches of a versioned
 * batch vector, we have to restore all previously removed single-pin and parallel nets
 * in order to uncontract the next batch vector. We create for each version of the
 * hypergraph a seperate batch uncontraction hierarchy (see createBatchUncontractionHierarchyOfVersion(...))
 */
// VersionedBatchVector DynamicGraph::createBatchUncontractionHierarchy(const size_t batch_size,
//                                                                           const bool test) {
// }

/**
 * Removes single-pin and parallel nets from the hypergraph. The weight
 * of a set of identical nets is aggregated in one representative hyperedge
 * and single-pin hyperedges are removed. Returns a vector of removed hyperedges.
 */
parallel::scalable_vector<ParallelHyperedge> DynamicGraph::removeSinglePinAndParallelHyperedges() {
}

/**
 * Restores a previously removed set of singple-pin and parallel hyperedges. Note, that hes_to_restore
 * must be exactly the same and given in the reverse order as returned by removeSinglePinAndParallelNets(...).
 */
void DynamicGraph::restoreSinglePinAndParallelNets(const parallel::scalable_vector<ParallelHyperedge>& hes_to_restore) {
  // Restores all previously removed hyperedges
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
// DynamicGraph::ContractionResult DynamicGraph::contract(const HypernodeID u,
//                                                                  const HypernodeID v,
//                                                                  const HypernodeWeight max_node_weight) {

//   // Acquire ownership in correct order to prevent deadlocks
//   if ( u < v ) {
//     acquireHypernode(u);
//     acquireHypernode(v);
//   } else {
//     acquireHypernode(v);
//     acquireHypernode(u);
//   }

//   // Contraction is valid if
//   //  1.) Contraction partner v is enabled
//   //  2.) There are no pending contractions on v
//   //  4.) Resulting node weight is less or equal than a predefined upper bound
//   const bool contraction_partner_valid =
//     nodeIsEnabled(v) && _contraction_tree.pendingContractions(v) == 0;
//   const bool less_or_equal_than_max_node_weight =
//     hypernode(u).weight() + hypernode(v).weight() <= max_node_weight;
//   if ( contraction_partner_valid && less_or_equal_than_max_node_weight ) {
//     ASSERT(nodeIsEnabled(u), "Hypernode" << u << "is disabled!");
//     hypernode(u).setWeight(nodeWeight(u) + nodeWeight(v));
//     hypernode(v).disable();
//     releaseHypernode(u);
//     releaseHypernode(v);

//     HypernodeID contraction_start = _contraction_index.load();
//     kahypar::ds::FastResetFlagArray<>& shared_incident_nets_u_and_v = _he_bitset.local();
//     shared_incident_nets_u_and_v.reset();
//     parallel::scalable_vector<HyperedgeID>& failed_hyperedge_contractions = _failed_hyperedge_contractions.local();
//     for ( const HyperedgeID& he : incidentEdges(v) ) {
//       // Try to acquire ownership of hyperedge. In case of success, we perform the
//       // contraction and otherwise, we remember the hyperedge and try later again.
//       if ( tryAcquireEdge(he) ) {
//         contractEdge(u, v, he, shared_incident_nets_u_and_v);
//         releaseEdge(he);
//       } else {
//         failed_hyperedge_contractions.push_back(he);
//       }
//     }

//     // Perform contraction on which we failed to acquire ownership on the first try
//     for ( const HyperedgeID& he : failed_hyperedge_contractions ) {
//       acquireEdge(he);
//       contractEdge(u, v, he, shared_incident_nets_u_and_v);
//       releaseEdge(he);
//     }

//     // Contract incident net lists of u and v
//     _incident_nets.contract(u, v, shared_incident_nets_u_and_v,
//       [&](const HypernodeID u) {
//         acquireHypernode(u);
//       }, [&](const HypernodeID u) {
//         releaseHypernode(u);
//       });
//     shared_incident_nets_u_and_v.reset();
//     failed_hyperedge_contractions.clear();

//     HypernodeID contraction_end = ++_contraction_index;
//     acquireHypernode(u);
//     _contraction_tree.unregisterContraction(u, v, contraction_start, contraction_end);
//     releaseHypernode(u);
//     return ContractionResult::CONTRACTED;
//   } else {
//     ContractionResult res = ContractionResult::PENDING_CONTRACTIONS;
//     if ( !less_or_equal_than_max_node_weight && nodeIsEnabled(v) &&
//          _contraction_tree.parent(v) == u ) {
//       _contraction_tree.unregisterContraction(u, v,
//         kInvalidHypernode, kInvalidHypernode, true /* failed */);
//       res = ContractionResult::WEIGHT_LIMIT_REACHED;
//     }
//     releaseHypernode(u);
//     releaseHypernode(v);
//     return res;
//   }
// }

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

/**
 * Computes a batch uncontraction hierarchy for a specific version of the hypergraph.
 * A batch is a vector of mementos (uncontractions) that are uncontracted in parallel.
 * Each time we perform single-pin and parallel net detection we create a new version of
 * the hypergraph.
 * A batch of uncontractions that is uncontracted in parallel must satisfy two conditions:
 *  1.) All representatives must be active vertices of the hypergraph
 *  2.) For a specific representative its contraction partners must be uncontracted in reverse
 *      contraction order. Meaning that a contraction (u, v) that happens before a contraction (u, w)
 *      must be uncontracted in a batch that is part of the same batch or a batch uncontracted after the
 *      batch which (u, w) is part of. This ensures that a parallel batch uncontraction does not
 *      increase the objective function.
 * We use the contraction tree to create a batch uncontraction order. Note, uncontractions from
 * different subtrees can be interleaved abitrary. To ensure condition 1.) we peform a BFS starting
 * from all roots of the contraction tree. Each BFS level induces a new batch. Since we contract
 * vertices in parallel its not possible to create a relative order of the contractions which is
 * neccassary for condition 2.). However, during a contraction we store a start and end "timepoint"
 * of a contraction. If two contractions time intervals do not intersect, we can determine
 * which contraction happens strictly before the other. If they intersect, it is not possible to
 * give a relative order. To ensure condition 2.) we sort the childs of a vertex in the contraction tree
 * after its time intervals. Once we add a uncontraction (u,v) to a batch, we also add all uncontractions
 * (u,w) to the batch which intersect the time interval of (u,v). To merge uncontractions of different
 * subtrees in a batch, we insert all eligble uncontractions into a max priority queue with the subtree
 * size of the contraction partner as key. We insert uncontractions into the current batch as long
 * as the maximum batch size is not reached or the PQ is empty. Once the batch reaches its maximum
 * batch size, we create a new empty batch. If the PQ is empty, we replace it with the PQ of the next
 * BFS level. With this approach heavy vertices are uncontracted earlier (subtree size in the PQ as key = weight of
 * a vertex for an unweighted hypergraph) such that average node weight of the hypergraph decreases faster and
 * local searches are more effective in early stages of the uncontraction hierarchy where hyperedge sizes are
 * usually smaller than on the original hypergraph.
 */

// BatchVector DynamicGraph::createBatchUncontractionHierarchyForVersion(BatchIndexAssigner& batch_assigner,
//                                                                            const size_t version) {

//   using PQ = std::priority_queue<PQBatchUncontractionElement,
//                                  parallel::scalable_vector<PQBatchUncontractionElement>,
//                                  PQElementComparator>;

//   // Checks if two contraction intervals intersect
//   auto does_interval_intersect = [&](const ContractionInterval& i1, const ContractionInterval& i2) {
//     if (i1.start == kInvalidHypernode || i2.start == kInvalidHypernode) {
//       return false;
//     }
//     return (i1.start <= i2.end && i1.end >= i2.end) ||
//             (i2.start <= i1.end && i2.end >= i1.end);
//   };

//   auto push_into_pq = [&](PQ& prio_q, const HypernodeID& u) {
//     auto it = _contraction_tree.childs(u);
//     auto current = it.begin();
//     auto end = it.end();
//     while ( current != end && _contraction_tree.version(*current) != version ) {
//       ++current;
//     }
//     if ( current != end ) {
//       prio_q.push(PQBatchUncontractionElement {
//         _contraction_tree.subtreeSize(*current), std::make_pair(current, end) } );
//     }
//   };

//   // Distribute roots of the contraction tree to local priority queues of
//   // each thread.
//   const size_t num_hardware_threads = std::thread::hardware_concurrency();
//   parallel::scalable_vector<PQ> local_pqs(num_hardware_threads);
//   const parallel::scalable_vector<HypernodeID>& roots = _contraction_tree.roots_of_version(version);
//   tbb::parallel_for(0UL, roots.size(), [&](const size_t i) {
//     const int cpu_id = sched_getcpu();
//     push_into_pq(local_pqs[cpu_id], roots[i]);
//   });

//   using LocalBatchAssignments = parallel::scalable_vector<BatchAssignment>;
//   parallel::scalable_vector<LocalBatchAssignments> local_batch_assignments(num_hardware_threads);
//   parallel::scalable_vector<size_t> local_batch_indices(num_hardware_threads, 0);
//   tbb::parallel_for(0UL, num_hardware_threads, [&](const size_t i) {
//     size_t& current_batch_index = local_batch_indices[i];
//     LocalBatchAssignments& batch_assignments = local_batch_assignments[i];
//     PQ& pq = local_pqs[i];
//     PQ next_pq;

//     while ( !pq.empty() ) {
//       // Iterator over the childs of a active vertex
//       auto it = pq.top()._iterator;
//       ASSERT(it.first != it.second);
//       const HypernodeID v = *it.first;
//       ASSERT(_contraction_tree.version(v) == version);
//       pq.pop();

//       const size_t start_idx = batch_assignments.size();
//       size_t num_uncontractions = 1;
//       const HypernodeID u = _contraction_tree.parent(v);
//       batch_assignments.push_back(BatchAssignment { u, v, 0UL, 0UL });
//       // Push contraction partner into pq for the next BFS level
//       push_into_pq(next_pq, v);

//       // Insert all childs of u that intersect the contraction time interval of
//       // (u,v) into the current batch
//       ++it.first;
//       ContractionInterval current_ival = _contraction_tree.interval(v);
//       while ( it.first != it.second && _contraction_tree.version(*it.first) == version ) {
//         const HypernodeID w = *it.first;
//         const ContractionInterval w_ival = _contraction_tree.interval(w);
//         if ( does_interval_intersect(current_ival, w_ival) ) {
//           ASSERT(_contraction_tree.parent(w) == u);
//           ++num_uncontractions;
//           batch_assignments.push_back(BatchAssignment { u, w, 0UL, 0UL });
//           current_ival.start = std::min(current_ival.start, w_ival.start);
//           current_ival.end = std::max(current_ival.end, w_ival.end);
//           push_into_pq(next_pq, w);
//         } else {
//           break;
//         }
//         ++it.first;
//       }

//       // If there are still childs left of u, we push the iterator again into the
//       // priority queue of the current BFS level.
//       if ( it.first != it.second && _contraction_tree.version(*it.first) == version ) {
//         pq.push(PQBatchUncontractionElement { _contraction_tree.subtreeSize(*it.first), it });
//       }

//       // Request batch index and its position within that batch
//       BatchAssignment assignment = batch_assigner.getBatchIndex(
//         current_batch_index, num_uncontractions);
//       for ( size_t j = start_idx; j < start_idx + num_uncontractions; ++j ) {
//         batch_assignments[j].batch_index = assignment.batch_index;
//         batch_assignments[j].batch_pos = assignment.batch_pos + (j - start_idx);
//       }
//       current_batch_index = assignment.batch_index;

//       if ( pq.empty() ) {
//         std::swap(pq, next_pq);
//         // Compute minimum batch index to which a thread assigned last.
//         // Afterwards, transmit information to batch assigner to speed up
//         // batch index computation.
//         ++current_batch_index;
//         size_t min_batch_index = current_batch_index;
//         for ( const size_t& batch_index : local_batch_indices ) {
//           min_batch_index = std::min(min_batch_index, batch_index);
//         }
//         batch_assigner.increaseHighWaterMark(min_batch_index);
//       }
//     }
//   });

//   ASSERT(verifyBatchIndexAssignments(batch_assigner, local_batch_assignments), "Batch asisignment failed");

//   // In the previous step we have calculated for each uncontraction a batch index and
//   // its position within that batch. We have to write the uncontractions
//   // into the global batch uncontraction vector.
//   const size_t num_batches = batch_assigner.numberOfNonEmptyBatches();
//   BatchVector batches(num_batches);
//   tbb::parallel_for(0UL, num_batches, [&](const size_t batch_index) {
//     batches[batch_index].resize(batch_assigner.batchSize(batch_index));
//   });

//   tbb::parallel_for(0UL, num_hardware_threads, [&](const size_t i) {
//     LocalBatchAssignments& batch_assignments = local_batch_assignments[i];
//     for ( const BatchAssignment& batch_assignment : batch_assignments ) {
//       const size_t batch_index = batch_assignment.batch_index;
//       const size_t batch_pos = batch_assignment.batch_pos;
//       ASSERT(batch_index < batches.size());
//       ASSERT(batch_pos < batches[batch_index].size());
//       batches[batch_index][batch_pos].u = batch_assignment.u;
//       batches[batch_index][batch_pos].v = batch_assignment.v;
//     }
//   });

//   while ( !batches.empty() && batches.back().empty() ) {
//     batches.pop_back();
//   }
//   std::reverse(batches.begin(), batches.end());

//   return batches;
// }

} // namespace ds
} // namespace mt_kahypar
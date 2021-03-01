/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/datastructures/dynamic_hypergraph.h"

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
void DynamicHypergraph::updateTotalWeight(const TaskGroupID) {
  _total_weight = tbb::parallel_reduce(tbb::blocked_range<HypernodeID>(ID(0), _num_hypernodes), 0,
    [this](const tbb::blocked_range<HypernodeID>& range, HypernodeWeight init) {
      HypernodeWeight weight = init;
      for (HypernodeID hn = range.begin(); hn < range.end(); ++hn) {
        if ( nodeIsEnabled(hn) ) {
          weight += this->_hypernodes[hn].weight();
        }
      }
      return weight;
    }, std::plus<HypernodeWeight>()) + _removed_degree_zero_hn_weight;
}

// ! Recomputes the total weight of the hypergraph (sequential)
void DynamicHypergraph::updateTotalWeight() {
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
bool DynamicHypergraph::registerContraction(const HypernodeID u, const HypernodeID v) {
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
size_t DynamicHypergraph::contract(const HypernodeID v,
                                   const HypernodeWeight max_node_weight) {
  ASSERT(_contraction_tree.parent(v) != v, "No contraction registered for hypernode" << v);

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
 * The two uncontraction functions are required by the partitioned hypergraph to restore
 * pin counts and gain cache values.
 */
void DynamicHypergraph::uncontract(const Batch& batch,
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

  _hes_to_resize_flag_array.reset();
  tbb::parallel_for(0UL, batch.size(), [&](const size_t i) {
    const Memento& memento = batch[i];
    ASSERT(!hypernode(memento.u).isDisabled(), "Hypernode" << memento.u << "is disabled");
    ASSERT(hypernode(memento.v).isDisabled(), "Hypernode" << memento.v << "is not invalid");

    // Restore incident net list of u and v
    const HypernodeID batch_index = hypernode(batch[0].v).batchIndex();
    _incident_nets.uncontract(memento.u, memento.v,
      [&](const HyperedgeID e) {
        // In that case, u and v were both previously part of hyperedge e.
        if ( !_hes_to_resize_flag_array[e] &&
             _hes_to_resize_flag_array.compare_and_set_to_true(e) ) {
          // This part is only triggered once for each hyperedge per batch uncontraction.
          // It restores all pins that are part of the current batch as contraction partners
          // in hyperedge e
          restoreHyperedgeSizeForBatch(e, batch_index);
        }

        acquireHyperedge(e);
        case_one_func(memento.u, memento.v, e);
        releaseHyperedge(e);
      }, [&](const HyperedgeID e) {
        // In that case only v was part of hyperedge e before and
        // u must be replaced by v in hyperedge e
        const size_t slot_of_u = findPositionOfPinInIncidenceArray(memento.u, e);

        acquireHyperedge(e);
        ASSERT(_incidence_array[slot_of_u] == memento.u);
        _incidence_array[slot_of_u] = memento.v;
        case_two_func(memento.u, memento.v, e);
        releaseHyperedge(e);
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
VersionedBatchVector DynamicHypergraph::createBatchUncontractionHierarchy(const size_t batch_size,
                                                                          const bool test) {
  const size_t num_versions = _version + 1;
  utils::Timer::instance().start_timer("finalize_contraction_tree", "Finalize Contraction Tree");
  // Finalizes the contraction tree such that it is traversable in a top-down fashion
  // and contains subtree size for each  tree node
  _contraction_tree.finalize(num_versions);
  utils::Timer::instance().stop_timer("finalize_contraction_tree");

  utils::Timer::instance().start_timer("create_versioned_batches", "Create Versioned Batches");
  VersionedBatchVector versioned_batches(num_versions);
  parallel::scalable_vector<size_t> batch_sizes_prefix_sum(num_versions, 0);
  BatchIndexAssigner batch_index_assigner(_num_hypernodes, batch_size);
  for ( size_t version = 0; version < num_versions; ++version ) {
    versioned_batches[version] =
      createBatchUncontractionHierarchyForVersion(
        batch_index_assigner, version);
    if ( version > 0 ) {
      batch_sizes_prefix_sum[version] =
        batch_sizes_prefix_sum[version - 1] + versioned_batches[version - 1].size();
    }
    batch_index_assigner.reset(versioned_batches[version].size());
  }
  utils::Timer::instance().stop_timer("create_versioned_batches");

  if ( !test ) {
    utils::Timer::instance().start_timer("prepare_hg_for_uncontraction", "Prepare HG For Uncontraction");

    // Store the batch index of each vertex in its hypernode data structure
    tbb::parallel_for(0UL, num_versions, [&](const size_t version) {
      tbb::parallel_for(0UL, versioned_batches[version].size(), [&](const size_t local_batch_idx) {
        const size_t batch_idx = batch_sizes_prefix_sum[version] + local_batch_idx;
        for ( const Memento& memento : versioned_batches[version][local_batch_idx] ) {
          hypernode(memento.v).setBatchIndex(batch_idx);
        }
      });
    });

    // Sort the invalid part of each hyperedge according to the batch indices of its pins
    tbb::parallel_for(ID(0), _num_hyperedges, [&](const HyperedgeID& he) {
      const size_t first_invalid_entry = hyperedge(he).firstInvalidEntry();
      const size_t last_invalid_entry = hyperedge(he + 1).firstEntry();
      std::sort(_incidence_array.begin() + first_invalid_entry,
                _incidence_array.begin() + last_invalid_entry,
                [&](const HypernodeID u, const HypernodeID v) {
                  ASSERT(hypernode(u).batchIndex() != std::numeric_limits<HypernodeID>::max(),
                    "Hypernode" << u << "is not contained in the uncontraction hierarchy");
                  ASSERT(hypernode(v).batchIndex() != std::numeric_limits<HypernodeID>::max(),
                    "Hypernode" << v << "is not contained in the uncontraction hierarchy");
                  return hypernode(u).batchIndex() > hypernode(v).batchIndex();
                });
    });
    utils::Timer::instance().stop_timer("prepare_hg_for_uncontraction");
  }

  return versioned_batches;
}

/**
 * Removes single-pin and parallel nets from the hypergraph. The weight
 * of a set of identical nets is aggregated in one representative hyperedge
 * and single-pin hyperedges are removed. Returns a vector of removed hyperedges.
 */
parallel::scalable_vector<ParallelHyperedge> DynamicHypergraph::removeSinglePinAndParallelHyperedges() {
  _removable_single_pin_and_parallel_nets.reset();
  // Remove singple-pin hyperedges directly from the hypergraph and
  // insert all other hyperedges into a bucket data structure such that
  // hyperedges with the same hash/footprint are placed in the same bucket.
  utils::Timer::instance().start_timer("preprocess_hyperedges", "Preprocess Hyperedges");
  StreamingVector<ParallelHyperedge> tmp_removed_hyperedges;
  ConcurrentBucketMap<ContractedHyperedgeInformation> hyperedge_hash_map;
  hyperedge_hash_map.reserve_for_estimated_number_of_insertions(_num_hyperedges);
  doParallelForAllEdges([&](const HyperedgeID& he) {
    const HypernodeID edge_size = edgeSize(he);
    if ( edge_size > 1 ) {
      const Hyperedge& e = hyperedge(he);
      const size_t footprint = e.hash();
      std::sort(_incidence_array.begin() + e.firstEntry(),
                _incidence_array.begin() + e.firstInvalidEntry());
      hyperedge_hash_map.insert(footprint,
        ContractedHyperedgeInformation { he, footprint, edge_size, true });
    } else {
      hyperedge(he).disable();
      _removable_single_pin_and_parallel_nets.set(he, true);
      tmp_removed_hyperedges.stream(ParallelHyperedge { he, kInvalidHyperedge });
    }
  });
  utils::Timer::instance().stop_timer("preprocess_hyperedges");

  // Helper function that checks if two hyperedges are parallel
  // Note, pins inside the hyperedges are sorted.
  auto check_if_hyperedges_are_parallel = [&](const HyperedgeID lhs,
                                              const HyperedgeID rhs) {
    const Hyperedge& lhs_he = hyperedge(lhs);
    const Hyperedge& rhs_he = hyperedge(rhs);
    if ( lhs_he.size() == rhs_he.size() ) {
      const size_t lhs_start = lhs_he.firstEntry();
      const size_t rhs_start = rhs_he.firstEntry();
      for ( size_t i = 0; i < lhs_he.size(); ++i ) {
        const size_t lhs_pos = lhs_start + i;
        const size_t rhs_pos = rhs_start + i;
        if ( _incidence_array[lhs_pos] != _incidence_array[rhs_pos] ) {
          return false;
        }
      }
      return true;
    } else {
      return false;
    }
  };

  // In the step before we placed hyperedges within a bucket data structure.
  // Hyperedges with the same hash/footprint are stored inside the same bucket.
  // We iterate now in parallel over each bucket and sort each bucket
  // after its hash. A bucket is processed by one thread and parallel
  // hyperedges are detected by comparing the pins of hyperedges with
  // the same hash.
  utils::Timer::instance().start_timer("parallel_net_detection", "Parallel Net Detection");
  tbb::parallel_for(0UL, hyperedge_hash_map.numBuckets(), [&](const size_t bucket) {
    auto& hyperedge_bucket = hyperedge_hash_map.getBucket(bucket);
    std::sort(hyperedge_bucket.begin(), hyperedge_bucket.end(),
      [&](const ContractedHyperedgeInformation& lhs, const ContractedHyperedgeInformation& rhs) {
        return lhs.hash < rhs.hash || (lhs.hash == rhs.hash && lhs.size < rhs.size)||
          (lhs.hash == rhs.hash && lhs.size == rhs.size && lhs.he < rhs.he);
      });

    // Parallel Hyperedge Detection
    for ( size_t i = 0; i < hyperedge_bucket.size(); ++i ) {
      ContractedHyperedgeInformation& contracted_he_lhs = hyperedge_bucket[i];
      if ( contracted_he_lhs.valid ) {
        const HyperedgeID lhs_he = contracted_he_lhs.he;
        HyperedgeWeight lhs_weight = hyperedge(lhs_he).weight();
        for ( size_t j = i + 1; j < hyperedge_bucket.size(); ++j ) {
          ContractedHyperedgeInformation& contracted_he_rhs = hyperedge_bucket[j];
          const HyperedgeID rhs_he = contracted_he_rhs.he;
          if ( contracted_he_rhs.valid &&
                contracted_he_lhs.hash == contracted_he_rhs.hash &&
                check_if_hyperedges_are_parallel(lhs_he, rhs_he) ) {
              // Hyperedges are parallel
              lhs_weight += hyperedge(rhs_he).weight();
              hyperedge(rhs_he).disable();
              _removable_single_pin_and_parallel_nets.set(rhs_he, true);
              contracted_he_rhs.valid = false;
              tmp_removed_hyperedges.stream( ParallelHyperedge { rhs_he, lhs_he } );
          } else if ( contracted_he_lhs.hash != contracted_he_rhs.hash  ) {
            // In case, hash of both are not equal we go to the next hyperedge
            // because we compared it with all hyperedges that had an equal hash
            break;
          }
        }
        hyperedge(lhs_he).setWeight(lhs_weight);
      }
    }
    hyperedge_hash_map.free(bucket);
  });
  utils::Timer::instance().stop_timer("parallel_net_detection");

  // Remove single-pin and parallel nets from incident net vector of vertices
  utils::Timer::instance().start_timer("postprocess_incident_nets", "Postprocess Incident Nets");
  doParallelForAllNodes([&](const HypernodeID& u) {
    _incident_nets.removeIncidentNets(u, _removable_single_pin_and_parallel_nets);
  });
  utils::Timer::instance().stop_timer("postprocess_incident_nets");

  utils::Timer::instance().start_timer("store_removed_hyperedges", "Store Removed Hyperedges");
  parallel::scalable_vector<ParallelHyperedge> removed_hyperedges = tmp_removed_hyperedges.copy_parallel();
  tmp_removed_hyperedges.clear_parallel();
  utils::Timer::instance().stop_timer("store_removed_hyperedges");

  ++_version;
  return removed_hyperedges;
}

/**
 * Restores a previously removed set of singple-pin and parallel hyperedges. Note, that hes_to_restore
 * must be exactly the same and given in the reverse order as returned by removeSinglePinAndParallelNets(...).
 */
void DynamicHypergraph::restoreSinglePinAndParallelNets(const parallel::scalable_vector<ParallelHyperedge>& hes_to_restore) {
  // Restores all previously removed hyperedges
  utils::Timer::instance().start_timer("restore_removed_nets", "Restore Removed Hyperedges");
  tbb::parallel_for(0UL, hes_to_restore.size(), [&](const size_t i) {
    const ParallelHyperedge& parallel_he = hes_to_restore[i];
    const HyperedgeID he = parallel_he.removed_hyperedge;
    ASSERT(!edgeIsEnabled(he), "Hyperedge" << he << "should be disabled");
    const bool is_parallel_net = parallel_he.representative != kInvalidHyperedge;
    hyperedge(he).enable();
    if ( is_parallel_net ) {
      const HyperedgeID rep = parallel_he.representative;
      ASSERT(edgeIsEnabled(rep), "Hyperedge" << rep << "should be enabled");
      Hyperedge& rep_he = hyperedge(rep);
      acquireHyperedge(rep);
      rep_he.setWeight(rep_he.weight() - hyperedge(he).weight());
      releaseHyperedge(rep);
    }
  });
  utils::Timer::instance().stop_timer("restore_removed_nets");

  utils::Timer::instance().start_timer("restore_in_incident_nets", "Restore HEs in Incident Nets");
  doParallelForAllNodes([&](const HypernodeID u) {
    _incident_nets.restoreIncidentNets(u);
  });
  utils::Timer::instance().stop_timer("restore_in_incident_nets");
  --_version;
}

// ! Copy dynamic hypergraph in parallel
DynamicHypergraph DynamicHypergraph::copy(const TaskGroupID task_group_id) {
  DynamicHypergraph hypergraph;

  hypergraph._num_hypernodes = _num_hypernodes;
  hypergraph._num_removed_hypernodes = _num_removed_hypernodes;
  hypergraph._removed_degree_zero_hn_weight = _removed_degree_zero_hn_weight;
  hypergraph._num_hyperedges = _num_hyperedges;
  hypergraph._num_removed_hyperedges = _num_removed_hyperedges;
  hypergraph._max_edge_size = _max_edge_size;
  hypergraph._num_pins = _num_pins;
  hypergraph._total_degree = _total_degree;
  hypergraph._total_weight = _total_weight;
  hypergraph._version = _version;
  hypergraph._contraction_index.store(_contraction_index.load());

  tbb::parallel_invoke([&] {
    hypergraph._hypernodes.resize(_hypernodes.size());
    memcpy(hypergraph._hypernodes.data(), _hypernodes.data(),
      sizeof(Hypernode) * _hypernodes.size());
  }, [&] {
    tbb::parallel_invoke([&] {
      hypergraph._incident_nets = _incident_nets.copy(task_group_id);
    }, [&] {
      hypergraph._acquired_hns.resize(_acquired_hns.size());
      tbb::parallel_for(ID(0), _num_hypernodes, [&](const HypernodeID& hn) {
        hypergraph._acquired_hns[hn] = _acquired_hns[hn];
      });
    });
  }, [&] {
    hypergraph._contraction_tree = _contraction_tree.copy(task_group_id);
  }, [&] {
    hypergraph._hyperedges.resize(_hyperedges.size());
    memcpy(hypergraph._hyperedges.data(), _hyperedges.data(),
      sizeof(Hyperedge) * _hyperedges.size());
  }, [&] {
    hypergraph._incidence_array.resize(_incidence_array.size());
    memcpy(hypergraph._incidence_array.data(), _incidence_array.data(),
      sizeof(HypernodeID) * _incidence_array.size());
  }, [&] {
    hypergraph._acquired_hes.resize(_num_hyperedges);
    tbb::parallel_for(ID(0), _num_hyperedges, [&](const HyperedgeID& he) {
      hypergraph._acquired_hes[he] = _acquired_hes[he];
    });
  }, [&] {
    hypergraph._hes_to_resize_flag_array =
      ThreadSafeFastResetFlagArray<>(_num_hyperedges);
  }, [&] {
    hypergraph._he_bitset = ThreadLocalBitset(_num_hyperedges);
  }, [&] {
    hypergraph._removable_single_pin_and_parallel_nets =
      kahypar::ds::FastResetFlagArray<>(_num_hyperedges);
  });
  return hypergraph;
}

// ! Copy dynamic hypergraph sequential
DynamicHypergraph DynamicHypergraph::copy() {
  DynamicHypergraph hypergraph;

  hypergraph._num_hypernodes = _num_hypernodes;
  hypergraph._num_removed_hypernodes = _num_removed_hypernodes;
  hypergraph._removed_degree_zero_hn_weight = _removed_degree_zero_hn_weight;
  hypergraph._num_hyperedges = _num_hyperedges;
  hypergraph._num_removed_hyperedges = _num_removed_hyperedges;
  hypergraph._max_edge_size = _max_edge_size;
  hypergraph._num_pins = _num_pins;
  hypergraph._total_degree = _total_degree;
  hypergraph._total_weight = _total_weight;
  hypergraph._version = _version;
  hypergraph._contraction_index.store(_contraction_index.load());

  hypergraph._hypernodes.resize(_hypernodes.size());
  memcpy(hypergraph._hypernodes.data(), _hypernodes.data(),
    sizeof(Hypernode) * _hypernodes.size());
  hypergraph._incident_nets = _incident_nets.copy();
  hypergraph._acquired_hns.resize(_num_hypernodes);
  for ( HypernodeID hn = 0; hn < _num_hypernodes; ++hn ) {
    hypergraph._acquired_hns[hn] = _acquired_hns[hn];
  }
  hypergraph._contraction_tree = _contraction_tree.copy();
  hypergraph._hyperedges.resize(_hyperedges.size());
  memcpy(hypergraph._hyperedges.data(), _hyperedges.data(),
    sizeof(Hyperedge) * _hyperedges.size());
  hypergraph._incidence_array.resize(_incidence_array.size());
  memcpy(hypergraph._incidence_array.data(), _incidence_array.data(),
    sizeof(HypernodeID) * _incidence_array.size());
  hypergraph._acquired_hes.resize(_num_hyperedges);
  for ( HyperedgeID he = 0; he < _num_hyperedges; ++he ) {
    hypergraph._acquired_hes[he] = _acquired_hes[he];
  }
  hypergraph._hes_to_resize_flag_array =
    ThreadSafeFastResetFlagArray<>(_num_hyperedges);
  hypergraph._he_bitset = ThreadLocalBitset(_num_hyperedges);
  hypergraph._removable_single_pin_and_parallel_nets =
    kahypar::ds::FastResetFlagArray<>(_num_hyperedges);

  return hypergraph;
}

void DynamicHypergraph::memoryConsumption(utils::MemoryTreeNode* parent) const {
  ASSERT(parent);

  parent->addChild("Hypernodes", sizeof(Hypernode) * _hypernodes.size());
  parent->addChild("Incident Nets", _incident_nets.size_in_bytes());
  parent->addChild("Hypernode Ownership Vector", sizeof(bool) * _acquired_hns.size());
  parent->addChild("Hyperedges", sizeof(Hyperedge) * _hyperedges.size());
  parent->addChild("Incidence Array", sizeof(HypernodeID) * _incidence_array.size());
  parent->addChild("Hyperedge Ownership Vector", sizeof(bool) * _acquired_hes.size());
  parent->addChild("Bitsets",
    ( _num_hyperedges * _he_bitset.size() ) / 8UL + sizeof(uint16_t) * _num_hyperedges);

  utils::MemoryTreeNode* contraction_tree_node = parent->addChild("Contraction Tree");
  _contraction_tree.memoryConsumption(contraction_tree_node);
}

// ! Only for testing
bool DynamicHypergraph::verifyIncidenceArrayAndIncidentNets() {
  bool success = true;
  tbb::parallel_invoke([&] {
    doParallelForAllNodes([&](const HypernodeID& hn) {
      for ( const HyperedgeID& he : incidentEdges(hn) ) {
        bool found = false;
        for ( const HypernodeID& pin : pins(he) ) {
          if ( pin == hn ) {
            found = true;
            break;
          }
        }
        if ( !found ) {
          LOG << "Hypernode" << hn << "not found in incidence array of net" << he;
          success = false;
        }
      }
    });
  }, [&] {
    doParallelForAllEdges([&](const HyperedgeID& he) {
      for ( const HypernodeID& pin : pins(he) ) {
        bool found = false;
        for ( const HyperedgeID& e : incidentEdges(pin) ) {
          if ( e == he ) {
            found = true;
            break;
          }
        }
        if ( !found ) {
          LOG << "Hyperedge" << he << "not found in incident nets of vertex" << pin;
          success = false;
        }
      }
    });
  });
  return success;
}

/**!
 * Contracts a previously registered contraction. The contraction of u and v is
 * performed if there are no pending contractions in the subtree of v and the
 * contractions respects the maximum allowed node weight. In case the contraction
 * was performed successfully, enum type CONTRACTED is returned. If contraction
 * was not performed either WEIGHT_LIMIT_REACHED (in case sum of both vertices is
 * greater than the maximum allowed node weight) or PENDING_CONTRACTIONS (in case
 * there are some unfinished contractions in the subtree of v) is returned.
 */
DynamicHypergraph::ContractionResult DynamicHypergraph::contract(const HypernodeID u,
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
    kahypar::ds::FastResetFlagArray<>& shared_incident_nets_u_and_v = _he_bitset.local();
    shared_incident_nets_u_and_v.reset();
    parallel::scalable_vector<HyperedgeID>& failed_hyperedge_contractions = _failed_hyperedge_contractions.local();
    for ( const HyperedgeID& he : incidentEdges(v) ) {
      // Try to acquire ownership of hyperedge. In case of success, we perform the
      // contraction and otherwise, we remember the hyperedge and try later again.
      if ( tryAcquireHyperedge(he) ) {
        contractHyperedge(u, v, he, shared_incident_nets_u_and_v);
        releaseHyperedge(he);
      } else {
        failed_hyperedge_contractions.push_back(he);
      }
    }

    // Perform contraction on which we failed to acquire ownership on the first try
    for ( const HyperedgeID& he : failed_hyperedge_contractions ) {
      acquireHyperedge(he);
      contractHyperedge(u, v, he, shared_incident_nets_u_and_v);
      releaseHyperedge(he);
    }

    // Contract incident net lists of u and v
    _incident_nets.contract(u, v, shared_incident_nets_u_and_v,
      [&](const HypernodeID u) {
        acquireHypernode(u);
      }, [&](const HypernodeID u) {
        releaseHypernode(u);
      });
    shared_incident_nets_u_and_v.reset();
    failed_hyperedge_contractions.clear();

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

// ! Performs the contraction of (u,v) inside hyperedge he
void DynamicHypergraph::contractHyperedge(const HypernodeID u,
                                          const HypernodeID v,
                                          const HyperedgeID he,
                                          kahypar::ds::FastResetFlagArray<>& shared_incident_nets_u_and_v) {
  Hyperedge& e = hyperedge(he);
  const HypernodeID pins_begin = e.firstEntry();
  const HypernodeID pins_end = e.firstInvalidEntry();
  HypernodeID slot_of_u = pins_end - 1;
  HypernodeID last_pin_slot = pins_end - 1;

  for (HypernodeID idx = pins_begin; idx != last_pin_slot; ++idx) {
    const HypernodeID pin = _incidence_array[idx];
    if (pin == v) {
      std::swap(_incidence_array[idx], _incidence_array[last_pin_slot]);
      --idx;
    } else if (pin == u) {
      slot_of_u = idx;
    }
  }

  ASSERT(_incidence_array[last_pin_slot] == v, "v is not last entry in adjacency array!");

  if (slot_of_u != last_pin_slot) {
    // Case 1:
    // Hyperedge e contains both u and v. Thus we don't need to connect u to e and
    // can just cut off the last entry in the edge array of e that now contains v.
    DBG << V(he) << ": Case 1";
    e.hash() -= kahypar::math::hash(v);
    e.decrementSize();
    shared_incident_nets_u_and_v.set(he, true);
  } else {
    DBG << V(he) << ": Case 2";
    // Case 2:
    // Hyperedge e does not contain u. Therefore we  have to connect e to the representative u.
    // This reuses the pin slot of v in e's incidence array (i.e. last_pin_slot!)
    e.hash() -= kahypar::math::hash(v);
    e.hash() += kahypar::math::hash(u);
    _incidence_array[last_pin_slot] = u;
  }
}

// ! Restore the size of the hyperedge to the size before the batch with
// ! index batch_index was contracted.
void DynamicHypergraph::restoreHyperedgeSizeForBatch(const HyperedgeID he,
                                                     const HypernodeID batch_index) {
  const size_t first_invalid_entry = hyperedge(he).firstInvalidEntry();
  const size_t last_invalid_entry = hyperedge(he + 1).firstEntry();
  const size_t edge_size = edgeSize(he);
  ASSERT(hypernode(_incidence_array[first_invalid_entry]).batchIndex() == batch_index);
  size_t pos = first_invalid_entry + 1;
  for ( ; pos < last_invalid_entry; ++pos ) {
    const HypernodeID pin = _incidence_array[pos];
    ASSERT(hypernode(pin).batchIndex() <= batch_index, V(he));
    if ( hypernode(pin).batchIndex() != batch_index ) {
      break;
    }
  }
  const size_t size_delta = pos - first_invalid_entry;
  hyperedge(he).setSize(edge_size + size_delta);
}

// ! Search for the position of pin u in hyperedge he in the incidence array
size_t DynamicHypergraph::findPositionOfPinInIncidenceArray(const HypernodeID u,
                                                            const HyperedgeID he) {
  const size_t first_valid_entry = hyperedge(he).firstEntry();
  const size_t first_invalid_entry = hyperedge(he).firstInvalidEntry();
  size_t slot_of_u = first_invalid_entry;
  for ( size_t pos = first_invalid_entry - 1; pos != first_valid_entry - 1; --pos ) {
    if ( u == _incidence_array[pos] ) {
      slot_of_u = pos;
      break;
    }
  }

  ASSERT(slot_of_u != first_invalid_entry,
    "Hypernode" << u << "is not incident to hyperedge" << he);
  return slot_of_u;
}

bool DynamicHypergraph::verifyBatchIndexAssignments(
  const BatchIndexAssigner& batch_assigner,
  const parallel::scalable_vector<parallel::scalable_vector<BatchAssignment>>& local_batch_assignments) const {
  parallel::scalable_vector<BatchAssignment> assignments;
  for ( size_t i = 0; i < local_batch_assignments.size(); ++i ) {
    for ( const BatchAssignment& batch_assign : local_batch_assignments[i] ) {
      assignments.push_back(batch_assign);
    }
  }
  std::sort(assignments.begin(), assignments.end(),
    [&](const BatchAssignment& lhs, const BatchAssignment& rhs) {
      return lhs.batch_index < rhs.batch_index ||
        (lhs.batch_index == rhs.batch_index && lhs.batch_pos < rhs.batch_pos);
    });

  if ( assignments.size() > 0 ) {
    if ( assignments[0].batch_index != 0 || assignments[0].batch_pos != 0 ) {
      LOG << "First uncontraction should start at batch 0 at position 0"
          << V(assignments[0].batch_index) << V(assignments[0].batch_pos);
      return false;
    }

    for ( size_t i = 1; i < assignments.size(); ++i ) {
      if ( assignments[i - 1].batch_index == assignments[i].batch_index ) {
        if ( assignments[i - 1].batch_pos + 1 != assignments[i].batch_pos ) {
          LOG << "Batch positions are not consecutive"
              << V(i) << V(assignments[i - 1].batch_pos) << V(assignments[i].batch_pos);
          return false;
        }
      } else {
        if ( assignments[i - 1].batch_index + 1 != assignments[i].batch_index ) {
          LOG << "Batch indices are not consecutive"
              << V(i) << V(assignments[i - 1].batch_index) << V(assignments[i].batch_index);
          return false;
        }
        if ( assignments[i].batch_pos != 0 ) {
          LOG << "First uncontraction of each batch should start at position 0"
              << V(assignments[i].batch_pos);
          return false;
        }
        if ( assignments[i - 1].batch_pos + 1 != batch_assigner.batchSize(assignments[i - 1].batch_index) ) {
          LOG << "Position of last uncontraction in batch" << assignments[i - 1].batch_index
              << "does not match size of batch"
              << V(assignments[i - 1].batch_pos) << V(batch_assigner.batchSize(assignments[i - 1].batch_index));
          return false;
        }
      }
    }
  }

  return true;
}

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

BatchVector DynamicHypergraph::createBatchUncontractionHierarchyForVersion(BatchIndexAssigner& batch_assigner,
                                                                           const size_t version) {

  using PQ = std::priority_queue<PQBatchUncontractionElement,
                                 parallel::scalable_vector<PQBatchUncontractionElement>,
                                 PQElementComparator>;

  // Checks if two contraction intervals intersect
  auto does_interval_intersect = [&](const ContractionInterval& i1, const ContractionInterval& i2) {
    if (i1.start == kInvalidHypernode || i2.start == kInvalidHypernode) {
      return false;
    }
    return (i1.start <= i2.end && i1.end >= i2.end) ||
            (i2.start <= i1.end && i2.end >= i1.end);
  };

  auto push_into_pq = [&](PQ& prio_q, const HypernodeID& u) {
    auto it = _contraction_tree.childs(u);
    auto current = it.begin();
    auto end = it.end();
    while ( current != end && _contraction_tree.version(*current) != version ) {
      ++current;
    }
    if ( current != end ) {
      prio_q.push(PQBatchUncontractionElement {
        _contraction_tree.subtreeSize(*current), std::make_pair(current, end) } );
    }
  };

  // Distribute roots of the contraction tree to local priority queues of
  // each thread.
  utils::Timer::instance().start_timer("assign_roots", "Assign Roots");
  const size_t num_hardware_threads = std::thread::hardware_concurrency();
  parallel::scalable_vector<PQ> local_pqs(num_hardware_threads);
  const parallel::scalable_vector<HypernodeID>& roots = _contraction_tree.roots_of_version(version);
  tbb::parallel_for(0UL, roots.size(), [&](const size_t i) {
    const int cpu_id = sched_getcpu();
    push_into_pq(local_pqs[cpu_id], roots[i]);
  });
  utils::Timer::instance().stop_timer("assign_roots");

  utils::Timer::instance().start_timer("create_tmp_batch_uncontraction_hierarchy", "Create Tmp Batch Unc. Hierarchy");
  using LocalBatchAssignments = parallel::scalable_vector<BatchAssignment>;
  parallel::scalable_vector<LocalBatchAssignments> local_batch_assignments(num_hardware_threads);
  parallel::scalable_vector<size_t> local_batch_indices(num_hardware_threads, 0);
  tbb::parallel_for(0UL, num_hardware_threads, [&](const size_t i) {
    size_t& current_batch_index = local_batch_indices[i];
    LocalBatchAssignments& batch_assignments = local_batch_assignments[i];
    PQ& pq = local_pqs[i];
    PQ next_pq;

    while ( !pq.empty() ) {
      // Iterator over the childs of a active vertex
      auto it = pq.top()._iterator;
      ASSERT(it.first != it.second);
      const HypernodeID v = *it.first;
      ASSERT(_contraction_tree.version(v) == version);
      pq.pop();

      const size_t start_idx = batch_assignments.size();
      size_t num_uncontractions = 1;
      const HypernodeID u = _contraction_tree.parent(v);
      batch_assignments.push_back(BatchAssignment { u, v, 0UL, 0UL });
      // Push contraction partner into pq for the next BFS level
      push_into_pq(next_pq, v);

      // Insert all childs of u that intersect the contraction time interval of
      // (u,v) into the current batch
      ++it.first;
      ContractionInterval current_ival = _contraction_tree.interval(v);
      while ( it.first != it.second && _contraction_tree.version(*it.first) == version ) {
        const HypernodeID w = *it.first;
        const ContractionInterval w_ival = _contraction_tree.interval(w);
        if ( does_interval_intersect(current_ival, w_ival) ) {
          ASSERT(_contraction_tree.parent(w) == u);
          ++num_uncontractions;
          batch_assignments.push_back(BatchAssignment { u, w, 0UL, 0UL });
          current_ival.start = std::min(current_ival.start, w_ival.start);
          current_ival.end = std::max(current_ival.end, w_ival.end);
          push_into_pq(next_pq, w);
        } else {
          break;
        }
        ++it.first;
      }

      // If there are still childs left of u, we push the iterator again into the
      // priority queue of the current BFS level.
      if ( it.first != it.second && _contraction_tree.version(*it.first) == version ) {
        pq.push(PQBatchUncontractionElement { _contraction_tree.subtreeSize(*it.first), it });
      }

      // Request batch index and its position within that batch
      BatchAssignment assignment = batch_assigner.getBatchIndex(
        current_batch_index, num_uncontractions);
      for ( size_t i = start_idx; i < start_idx + num_uncontractions; ++i ) {
        batch_assignments[i].batch_index = assignment.batch_index;
        batch_assignments[i].batch_pos = assignment.batch_pos + (i - start_idx);
      }
      current_batch_index = assignment.batch_index;

      if ( pq.empty() ) {
        std::swap(pq, next_pq);
        // Compute minimum batch index to which a thread assigned last.
        // Afterwards, transmit information to batch assigner to speed up
        // batch index computation.
        ++current_batch_index;
        size_t min_batch_index = current_batch_index;
        for ( const size_t& batch_index : local_batch_indices ) {
          min_batch_index = std::min(min_batch_index, batch_index);
        }
        batch_assigner.increaseHighWaterMark(min_batch_index);
      }
    }
  });
  utils::Timer::instance().stop_timer("create_tmp_batch_uncontraction_hierarchy");

  ASSERT(verifyBatchIndexAssignments(batch_assigner, local_batch_assignments), "Batch asisignment failed");

  // In the previous step we have calculated for each uncontraction a batch index and
  // its position within that batch. We have to write the uncontractions
  // into the global batch uncontraction vector.
  utils::Timer::instance().start_timer("create_batch_vector", "Create Batch Vector");
  const size_t num_batches = batch_assigner.numberOfNonEmptyBatches();
  BatchVector batches(num_batches);
  tbb::parallel_for(0UL, num_batches, [&](const size_t batch_index) {
    batches[batch_index].resize(batch_assigner.batchSize(batch_index));
  });

  tbb::parallel_for(0UL, num_hardware_threads, [&](const size_t i) {
    LocalBatchAssignments& batch_assignments = local_batch_assignments[i];
    for ( const BatchAssignment& batch_assignment : batch_assignments ) {
      const size_t batch_index = batch_assignment.batch_index;
      const size_t batch_pos = batch_assignment.batch_pos;
      ASSERT(batch_index < batches.size());
      ASSERT(batch_pos < batches[batch_index].size());
      batches[batch_index][batch_pos].u = batch_assignment.u;
      batches[batch_index][batch_pos].v = batch_assignment.v;
    }
  });
  utils::Timer::instance().stop_timer("create_batch_vector");

  while ( !batches.empty() && batches.back().empty() ) {
    batches.pop_back();
  }
  std::reverse(batches.begin(), batches.end());

  return batches;
}

} // namespace ds
} // namespace mt_kahypar
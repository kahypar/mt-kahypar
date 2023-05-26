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

#include "mt-kahypar/partition/refinement/gains/process_mapping/process_mapping_gain_cache.h"

#include "tbb/parallel_for.h"
#include "tbb/enumerable_thread_specific.h"
#include "tbb/concurrent_vector.h"

#include "mt-kahypar/definitions.h"

namespace mt_kahypar {

template<typename PartitionedHypergraph>
void ProcessMappingGainCache::initializeGainCache(const PartitionedHypergraph& partitioned_hg) {
  ASSERT(!_is_initialized, "Gain cache is already initialized");
  ASSERT(_k == kInvalidPartition || _k == partitioned_hg.k(), "Gain cache was already initialized for a different k");
  allocateGainTable(partitioned_hg.topLevelNumNodes(), partitioned_hg.topLevelNumEdges(), partitioned_hg.k());
  initializeAdjacentBlocksOfEachNode(partitioned_hg);

  // Compute gain of all nodes
  tbb::parallel_for(tbb::blocked_range<HypernodeID>(HypernodeID(0), partitioned_hg.initialNumNodes()),
    [&](tbb::blocked_range<HypernodeID>& r) {
      vec<HyperedgeWeight>& benefit_aggregator = _ets_benefit_aggregator.local();
      for (HypernodeID u = r.begin(); u < r.end(); ++u) {
        if ( partitioned_hg.nodeIsEnabled(u)) {
          initializeGainCacheEntryForNode(partitioned_hg, u, benefit_aggregator);
        }
      }
    });

  _is_initialized = true;
}

bool ProcessMappingGainCache::triggersDeltaGainUpdate(const SyncronizedEdgeUpdate& sync_update) {
  return sync_update.pin_count_in_from_part_after == 0 ||
         sync_update.pin_count_in_from_part_after == 1 ||
         sync_update.pin_count_in_to_part_after == 1 ||
         sync_update.pin_count_in_to_part_after == 2;
}

void ProcessMappingGainCache::updateVersionOfHyperedge(const SyncronizedEdgeUpdate& sync_update) {
  if ( triggersDeltaGainUpdate(sync_update) ) {
    ASSERT(UL(sync_update.he) < _version.size());
    // The move will induce a gain cache update. In this case, we increment the version ID
    // of the hyperedge such that concurrent running initializations of gain entries are
    // notified and rerun the initialization step if this affected the gain computation.
    ++_version[sync_update.he].version;
  }
}

namespace {
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
HyperedgeWeight gainOfHyperedge(const PartitionID from,
                                const PartitionID to,
                                const HyperedgeWeight edge_weight,
                                const ProcessGraph& process_graph,
                                ds::PinCountSnapshot& pin_counts,
                                ds::Bitset& connectivity_set) {
  const HypernodeID pin_count_in_from_part = pin_counts.pinCountInPart(from);
  const HyperedgeWeight current_distance = process_graph.distance(connectivity_set);
  if ( pin_count_in_from_part == 1 ) {
    ASSERT(connectivity_set.isSet(from));
    connectivity_set.unset(from);
  }
  const HyperedgeWeight distance_with_to = process_graph.distanceWithBlock(connectivity_set, to);
  if ( pin_count_in_from_part == 1 ) {
    connectivity_set.set(from);
  }
  return (current_distance - distance_with_to) * edge_weight;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
void reconstructConnectivitySetAndPinCountsBeforeMove(const SyncronizedEdgeUpdate& sync_update,
                                                      ds::Bitset& connectivity_set,
                                                      ds::PinCountSnapshot& pin_counts) {
  if ( sync_update.pin_count_in_from_part_after == 0 ) {
    ASSERT(!connectivity_set.isSet(sync_update.from));
    connectivity_set.set(sync_update.from);
  }
  if ( sync_update.pin_count_in_to_part_after == 1 ) {
    ASSERT(connectivity_set.isSet(sync_update.to));
    connectivity_set.unset(sync_update.to);
  }
  pin_counts.setPinCountInPart(sync_update.from, sync_update.pin_count_in_from_part_after + 1);
  pin_counts.setPinCountInPart(sync_update.to, sync_update.pin_count_in_to_part_after - 1);
}
}

template<typename PartitionedHypergraph>
void ProcessMappingGainCache::deltaGainUpdate(const PartitionedHypergraph& partitioned_hg,
                                              const SyncronizedEdgeUpdate& sync_update) {
  ASSERT(_is_initialized, "Gain cache is not initialized");
  ASSERT(sync_update.connectivity_set_after);
  ASSERT(sync_update.pin_counts_after);
  ASSERT(sync_update.process_graph);

  if ( triggersDeltaGainUpdate(sync_update) ) {
    const HyperedgeID he = sync_update.he;
    const PartitionID from = sync_update.from;
    const PartitionID to = sync_update.to;
    const HyperedgeWeight edge_weight = sync_update.edge_weight;
    const HypernodeID pin_count_in_from_part_after = sync_update.pin_count_in_from_part_after;
    const HypernodeID pin_count_in_to_part_after = sync_update.pin_count_in_to_part_after;
    const ProcessGraph& process_graph = *sync_update.process_graph;
    ds::Bitset& connectivity_set = *sync_update.connectivity_set_after;
    ds::PinCountSnapshot& pin_counts = *sync_update.pin_counts_after;

    if ( pin_count_in_from_part_after == 0 || pin_count_in_to_part_after == 1 ) {
      // Connectivity set has changed
      // => Recompute gain of hyperedge for all pins and their adjacent blocks

      // Compute new gain of hyperedge for all pins and their adjacent blocks and
      // add it to the gain cache entries
      for ( const HypernodeID& pin : partitioned_hg.pins(he) ) {
        const PartitionID source = partitioned_hg.partID(pin);
        for ( const PartitionID& target : _adjacent_blocks.connectivitySet(pin) ) {
          const HyperedgeWeight gain_after = gainOfHyperedge(
            source, target, edge_weight, process_graph, pin_counts, connectivity_set);
          _gain_cache[benefit_index(pin, target)].add_fetch(gain_after, std::memory_order_relaxed);
        }
      }

      // Reconstruct connectivity set and pin counts before the node move
      reconstructConnectivitySetAndPinCountsBeforeMove(sync_update, connectivity_set, pin_counts);
      // Compute old gain of hyperedge for all pins and their adjacent blocks and
      // subtract it from the gain cache entries
      for ( const HypernodeID& pin : partitioned_hg.pins(he) ) {
        const PartitionID source = partitioned_hg.partID(pin);
        for ( const PartitionID& target : _adjacent_blocks.connectivitySet(pin) ) {
          const HyperedgeWeight gain_before = gainOfHyperedge(
            source, target, edge_weight, process_graph, pin_counts, connectivity_set);
          _gain_cache[benefit_index(pin, target)].sub_fetch(gain_before, std::memory_order_relaxed);
        }
      }
    } else {
      if ( pin_count_in_from_part_after == 1 ) {
        // In this case, there is only one pin left in block `from` and moving it to another block
        // would remove the block from the connectivity set. Thus, we search for the last remaining pin
        // in that block and update its gains for moving it to all its adjacent blocks.
        for ( const HypernodeID& u : partitioned_hg.pins(he) ) {
          if ( partitioned_hg.partID(u) == from ) {
            for ( const PartitionID& target : _adjacent_blocks.connectivitySet(u) ) {
              if ( from != target ) {
                // Compute new gain of hyperedge for moving u to the target block
                const HyperedgeWeight gain = gainOfHyperedge(
                  from, target, edge_weight, process_graph, pin_counts, connectivity_set);
                _gain_cache[benefit_index(u, target)].add_fetch(gain, std::memory_order_relaxed);

                // Before the node move, we would have increase the connectivity of the hyperedge
                // if we would have moved u to a block not in the connectivity set of the hyperedge.
                // Thus, we subtract the old gain from gain cache entry.
                const HypernodeID pin_count_target_part_before = target == to ?
                  pin_count_in_to_part_after - 1 : pin_counts.pinCountInPart(target);
                if ( pin_count_target_part_before == 0 ) {
                  // The target part was not part of the connectivity set of the hyperedge before the move.
                  // Thus, moving u to that block would have increased the connectivity of the hyperedge.
                  // However, this is no longer the case since moving u out of its block would remove the
                  // block from the connectivity set.
                  const bool was_set = connectivity_set.isSet(target);
                  connectivity_set.unset(target);
                  const HyperedgeWeight distance_before = process_graph.distance(connectivity_set);
                  const HyperedgeWeight distance_after = process_graph.distanceWithBlock(connectivity_set, target);
                  const HyperedgeWeight gain_before = (distance_before - distance_after) * edge_weight;
                  _gain_cache[benefit_index(u, target)].sub_fetch(gain_before, std::memory_order_relaxed);
                  if ( was_set ) connectivity_set.set(target);
                }
              }
            }
          }
        }
      }

      if (pin_count_in_to_part_after == 2) {
        // In this case, there are now two pins in block `to`. However, moving out the previously last pin
        // of block `to` would have decreased the connectivity of the hyperedge. This is no longer the case
        // since there are two pins in the block. Thus, we search for this pin and update its gain.
        for ( const HypernodeID& u : partitioned_hg.pins(he) ) {
          if ( partitioned_hg.partID(u) == to ) {
            for ( const PartitionID& target : _adjacent_blocks.connectivitySet(u) ) {
              if ( target != to ) {
                // Compute new gain of hyperedge for moving u to the target block
                const HyperedgeWeight gain = gainOfHyperedge(
                  to, target, edge_weight, process_graph, pin_counts, connectivity_set);
                _gain_cache[benefit_index(u, target)].add_fetch(gain, std::memory_order_relaxed);

                // Before the node move, we would have decreased the connectivity of the hyperedge
                // if we would have moved u to a block in the connecivity set or replaced its block
                // with another if we would have moved it to block not in the connectivity set.
                // Thus, we subtract the old gain from gain cache entry.
                const HypernodeID pin_count_target_part_before = target == from ?
                  pin_count_in_from_part_after + 1 : pin_counts.pinCountInPart(target);
                const bool was_set = connectivity_set.isSet(target);
                if ( pin_count_target_part_before == 0 ) connectivity_set.unset(target);
                const HyperedgeWeight distance_before = process_graph.distance(connectivity_set);
                HyperedgeWeight distance_after = 0;
                if ( pin_count_target_part_before > 0 ) {
                  // The target block was part of the connectivity set before the node move.
                  // Thus, moving u out of its block would have decreased the connectivity of
                  // the hyperedge.
                  distance_after = process_graph.distanceWithoutBlock(connectivity_set, to);
                } else {
                  // The target block was not part of the connectivity set before the node move.
                  // Thus, moving u out of its block would have replaced block `to` with the target block
                  // in the connectivity set.
                  distance_after = process_graph.distanceAfterExchangingBlocks(connectivity_set, to, target);
                }
                const HyperedgeWeight gain_before = (distance_before - distance_after) * edge_weight;
                _gain_cache[benefit_index(u, target)].sub_fetch(gain_before, std::memory_order_relaxed);
                if ( was_set ) connectivity_set.set(target);
              }
            }
          }
        }
      }
    }

    // Update gain version of hyperedge. If the update version is equal to the version
    // of the hyperedge, then we know that all gain cache updates are completed. This is
    // important for initializing gain entries while simultanously running gain cache updates.
    ++_version[sync_update.he].update_version;
  }

  // We update the adjacent blocks of nodes affected by this gain cache update
  // which will then trigger initialization of gain entries if a node becomes adjacent
  // to a new block.
  updateAdjacentBlocks(partitioned_hg, sync_update);
}

template<typename PartitionedHypergraph>
void ProcessMappingGainCache::uncontractUpdateAfterRestore(const PartitionedHypergraph&,
                                                           const HypernodeID,
                                                           const HypernodeID,
                                                           const HyperedgeID,
                                                           const HypernodeID) {
  if ( _is_initialized ) {

  }
}

template<typename PartitionedHypergraph>
void ProcessMappingGainCache::uncontractUpdateAfterReplacement(const PartitionedHypergraph&,
                                                               const HypernodeID,
                                                               const HypernodeID,
                                                               const HyperedgeID) {
  // In this case, u is replaced by v in hyperedge he
  // => Pin counts of hyperedge he does not change
  if ( _is_initialized ) {

  }
}

void ProcessMappingGainCache::restoreSinglePinHyperedge(const HypernodeID,
                                                        const PartitionID,
                                                        const HyperedgeWeight) {
  if ( _is_initialized ) {

  }
}

template<typename PartitionedHypergraph>
void ProcessMappingGainCache::initializeAdjacentBlocksOfEachNode(const PartitionedHypergraph& partitioned_hg) {
  // Initialize adjacent blocks of each node
  partitioned_hg.doParallelForAllNodes([&](const HypernodeID& hn) {
    _adjacent_blocks.clear(hn);
    for ( PartitionID to = 0; to < _k; ++to ) {
      _num_incident_edges_of_block[benefit_index(hn, to)].store(0, std::memory_order_relaxed);
    }
    for ( const HyperedgeID& he : partitioned_hg.incidentEdges(hn) ) {
      for ( const PartitionID& block : partitioned_hg.connectivitySet(he) ) {
        incrementIncidentEdges(hn, block);
      }
    }
  });
}

template<typename PartitionedHypergraph>
void ProcessMappingGainCache::updateAdjacentBlocks(const PartitionedHypergraph& partitioned_hg,
                                                   const SyncronizedEdgeUpdate& sync_update) {
  if ( sync_update.pin_count_in_from_part_after == 0 ) {
    // The node move has removed the source block of the move from the
    // connectivity set of the hyperedge. We therefore decrement the number of
    // incident edges in the source block for each pin of the hyperedge. If this
    // decreases the counter to zero for some pin, we remove the source block
    // from the adjacent blocks of that pin.
    for ( const HypernodeID& pin : partitioned_hg.pins(sync_update.he) ) {
      decrementIncidentEdges(pin, sync_update.from);
    }
  }
  if ( sync_update.pin_count_in_to_part_after == 1 ) {
    // The node move has added the target block of the move to the
    // connectivity set of the hyperedge. We therefore increment the number of
    // incident edges in the target block for each pin of the hyperedge. If this
    // increases the counter to one for some pin, we add the target block
    // to the adjacent blocks of that pin. Moreover, since we only compute gain
    // cache entries to adjacent blocks, we initialize the gain cache entry
    // for that pin and target block.
    for ( const HypernodeID& pin : partitioned_hg.pins(sync_update.he) ) {
      const HyperedgeID incident_edges_after = incrementIncidentEdges(pin, sync_update.to);
      if ( incident_edges_after == 1 ) {
        ASSERT(sync_update.edge_locks);
        initializeGainCacheEntry(partitioned_hg, pin, sync_update.to, *sync_update.edge_locks);
      }
    }
  }
}

HyperedgeID ProcessMappingGainCache::incrementIncidentEdges(const HypernodeID u, const PartitionID to) {
  const HyperedgeID incident_count_after =
    _num_incident_edges_of_block[benefit_index(u, to)].add_fetch(1, std::memory_order_relaxed);
  if ( incident_count_after == 1 ) {
    ASSERT(!_adjacent_blocks.contains(u, to));
    _gain_cache[benefit_index(u, to)].store(0, std::memory_order_relaxed);
    _adjacent_blocks.add(u, to);
  }
  return incident_count_after;
}

HyperedgeID ProcessMappingGainCache::decrementIncidentEdges(const HypernodeID u, const PartitionID to) {
  const HyperedgeID incident_count_after =
    _num_incident_edges_of_block[benefit_index(u, to)].sub_fetch(1, std::memory_order_relaxed);
  if ( incident_count_after == 0 ) {
    ASSERT(_adjacent_blocks.contains(u, to));
    _adjacent_blocks.remove(u, to);
  }
  return incident_count_after;
}

template<typename PartitionedHypergraph>
void ProcessMappingGainCache::initializeGainCacheEntryForNode(const PartitionedHypergraph& partitioned_hg,
                                                              const HypernodeID u,
                                                              vec<Gain>& benefit_aggregator) {
  ASSERT(partitioned_hg.hasProcessGraph());
  const ProcessGraph& process_graph = *partitioned_hg.processGraph();
  const PartitionID from = partitioned_hg.partID(u);

  // We only compute the gain to adjacent blocks of a node and initialize them here.
  // The gain to non-adjacent blocks is -inf.
  for ( const PartitionID& to : _adjacent_blocks.connectivitySet(u) ) {
    benefit_aggregator[to] = 0;
  }

  for ( const HyperedgeID& he : partitioned_hg.incidentEdges(u) ) {
    const HyperedgeWeight edge_weight = partitioned_hg.edgeWeight(he);
    ds::Bitset& connectivity_set = partitioned_hg.deepCopyOfConnectivitySet(he);

    const HyperedgeWeight current_distance = process_graph.distance(connectivity_set);
    if ( partitioned_hg.pinCountInPart(he, from) == 1 ) {
      // Moving the node out of its current block removes
      // its block from the connectivity set
      connectivity_set.unset(from);
    }
    // Compute gain to all adjacent blocks
    for ( const PartitionID& to : _adjacent_blocks.connectivitySet(u) ) {
      const HyperedgeWeight distance_with_to =
        process_graph.distanceWithBlock(connectivity_set, to);
      benefit_aggregator[to] += ( current_distance - distance_with_to ) * edge_weight;
    }
  }

  for ( PartitionID to = 0; to < _k; ++to ) {
    _gain_cache[benefit_index(u, to)].store(benefit_aggregator[to], std::memory_order_relaxed);
    benefit_aggregator[to] = std::numeric_limits<Gain>::min();
  }
}



template<typename PartitionedHypergraph>
void ProcessMappingGainCache::initializeGainCacheEntry(const PartitionedHypergraph& partitioned_hg,
                                                       const HypernodeID hn,
                                                       const PartitionID to,
                                                       ds::Array<SpinLock>& edge_locks) {
  ASSERT(partitioned_hg.hasProcessGraph());
  const ProcessGraph& process_graph = *partitioned_hg.processGraph();
  const HypernodeID from = partitioned_hg.partID(hn);
  vec<uint32_t>& seen_versions = _ets_version.local();
  bool success = false;
  while ( !success ) {
    success = true;
    seen_versions.clear();
    HyperedgeWeight gain = 0;
    for ( const HyperedgeID& he : partitioned_hg.incidentEdges(hn) ) {
      edge_locks[partitioned_hg.uniqueEdgeID(he)].lock();
      // The internal data structures in the partitioned hypergraph are updated
      // in one transaction and each update is assoicated with a version ID. We
      // retrieve here the actual state of the connectivity set of the hyperedge
      // with its version ID. If this version ID changes after the gain computation,
      // we know that we computed the gain on outdated information and retry.
      const uint32_t he_version = _version[he].version.load(std::memory_order_relaxed);
      ds::Bitset& connectivity_set = partitioned_hg.deepCopyOfConnectivitySet(he);
      edge_locks[partitioned_hg.uniqueEdgeID(he)].unlock();

      const uint32_t update_version = _version[he].update_version.load(std::memory_order_relaxed);
      ASSERT(update_version <= he_version);
      if ( update_version < he_version ) {
        // There are still pending gain cache updates that must be finished
        // before we initialize the gain cache entry.
        success = false;
        break;
      }
      seen_versions.push_back(he_version);

      // Now compute gain of moving node hn to block `to` for hyperedge
      const HyperedgeWeight current_distance = process_graph.distance(connectivity_set);
      if ( partitioned_hg.pinCountInPart(he, from) == 1 ) {
        // Moving the node out of its current block removes
        // its block from the connectivity set
        connectivity_set.unset(from);
      }
      const HyperedgeWeight distance_with_to =
        process_graph.distanceWithBlock(connectivity_set, to);
      gain += (current_distance - distance_with_to) * partitioned_hg.edgeWeight(he);
    }
    _gain_cache[benefit_index(hn, to)].store(gain, std::memory_order_relaxed);

    // Check if versions of an incident hyperedge has changed in the meantime.
    // If not, gain cache entry is correct. Otherwise, recompute it.
    if ( success ) {
      ASSERT(seen_versions.size() == UL(partitioned_hg.nodeDegree(hn)),
        V(hn) << V(seen_versions.size()) << V(partitioned_hg.nodeDegree(hn)));
      size_t idx = 0;
      for ( const HyperedgeID& he : partitioned_hg.incidentEdges(hn) ) {
        if ( seen_versions[idx++] != _version[he].version.load(std::memory_order_relaxed) ) {
          success = false;
          break;
        }
      }
    }
  }
}

namespace {
#define PROCESS_MAPPING_INITIALIZE_GAIN_CACHE(X) void ProcessMappingGainCache::initializeGainCache(const X&)
#define PROCESS_MAPPING_DELTA_GAIN_UPDATE(X) void ProcessMappingGainCache::deltaGainUpdate(const X&,                     \
                                                                                           const SyncronizedEdgeUpdate&)
#define PROCESS_MAPPING_RESTORE_UPDATE(X) void ProcessMappingGainCache::uncontractUpdateAfterRestore(const X&,          \
                                                                                                     const HypernodeID, \
                                                                                                     const HypernodeID, \
                                                                                                     const HyperedgeID, \
                                                                                                     const HypernodeID)
#define PROCESS_MAPPING_REPLACEMENT_UPDATE(X) void ProcessMappingGainCache::uncontractUpdateAfterReplacement(const X&,            \
                                                                                                             const HypernodeID,   \
                                                                                                             const HypernodeID,   \
                                                                                                             const HyperedgeID)
#define PROCESS_MAPPING_INIT_ADJACENT_BLOCKS(X) void ProcessMappingGainCache::initializeAdjacentBlocksOfEachNode(const X&)
#define PROCESS_MAPPING_UPDATE_ADJACENT_BLOCKS(X) void ProcessMappingGainCache::updateAdjacentBlocks(const X&,                     \
                                                                                                     const SyncronizedEdgeUpdate&)
#define PROCESS_MAPPING_INIT_GAIN_CACHE_ENTRY(X) void ProcessMappingGainCache::initializeGainCacheEntryForNode(const X&,           \
                                                                                                               const HypernodeID,  \
                                                                                                               vec<Gain>&)
#define PROCESS_MAPPING_INIT_LAZY_GAIN_CACHE_ENTRY(X) void ProcessMappingGainCache::initializeGainCacheEntry(const X&,             \
                                                                                                             const HypernodeID,    \
                                                                                                             const PartitionID,    \
                                                                                                             ds::Array<SpinLock>&)
}

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(PROCESS_MAPPING_INITIALIZE_GAIN_CACHE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(PROCESS_MAPPING_DELTA_GAIN_UPDATE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(PROCESS_MAPPING_RESTORE_UPDATE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(PROCESS_MAPPING_REPLACEMENT_UPDATE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(PROCESS_MAPPING_INIT_ADJACENT_BLOCKS)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(PROCESS_MAPPING_UPDATE_ADJACENT_BLOCKS)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(PROCESS_MAPPING_INIT_GAIN_CACHE_ENTRY)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(PROCESS_MAPPING_INIT_LAZY_GAIN_CACHE_ENTRY)

}  // namespace mt_kahypar

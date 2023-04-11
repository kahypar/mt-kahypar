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

#include "kahypar/meta/typelist.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/partition/refinement/gains/km1/km1_gain_cache.h"
#include "mt-kahypar/partition/refinement/gains/cut/cut_gain_cache.h"
#ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
#include "mt-kahypar/partition/refinement/gains/cut_for_graphs/cut_gain_cache_for_graphs.h"
#endif
#include "mt-kahypar/macros.h"

namespace mt_kahypar {

struct gain_cache_s;
typedef struct  {
  gain_cache_s* gain_cache;
  GainPolicy type;
} gain_cache_t;

class GainCachePtr {

 public:
  static gain_cache_t constructGainCache(const GainPolicy& type) {
    switch(type) {
      case GainPolicy::km1: return constructGainCache<Km1GainCache>();
      case GainPolicy::cut: return constructGainCache<CutGainCache>();
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case GainPolicy::cut_for_graphs: return constructGainCache<GraphCutGainCache>();
      #endif
      case GainPolicy::none:
        ERR("No gain policy set");
    }
    return gain_cache_t { nullptr, GainPolicy::none };
  }

  static void deleteGainCache(gain_cache_t gain_cache) {
    if ( gain_cache.gain_cache ) {
      switch(gain_cache.type) {
        case GainPolicy::cut:
          delete reinterpret_cast<CutGainCache*>(gain_cache.gain_cache); break;
        case GainPolicy::km1:
          delete reinterpret_cast<Km1GainCache*>(gain_cache.gain_cache); break;
        #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
        case GainPolicy::cut_for_graphs:
          delete reinterpret_cast<GraphCutGainCache*>(gain_cache.gain_cache); break;
        #endif
        case GainPolicy::none: break;
      }
    }
  }

  template<typename PartitionedHypergraph>
  static void initializeGainCache(const PartitionedHypergraph& partitioned_hg,
                                  gain_cache_t gain_cache) {
    switch(gain_cache.type) {
      case GainPolicy::cut:
        cast<CutGainCache>(gain_cache).initializeGainCache(partitioned_hg); break;
      case GainPolicy::km1:
        cast<Km1GainCache>(gain_cache).initializeGainCache(partitioned_hg); break;
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case GainPolicy::cut_for_graphs:
        cast<GraphCutGainCache>(gain_cache).initializeGainCache(partitioned_hg); break;
      #endif
      case GainPolicy::none: break;
    }
  }

  static void resetGainCache(gain_cache_t gain_cache) {
    switch(gain_cache.type) {
      case GainPolicy::cut:
        cast<CutGainCache>(gain_cache).reset(); break;
      case GainPolicy::km1:
        cast<Km1GainCache>(gain_cache).reset(); break;
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case GainPolicy::cut_for_graphs:
        cast<GraphCutGainCache>(gain_cache).reset(); break;
      #endif
      case GainPolicy::none: break;
    }
  }

  template<typename PartitionedHypergraph>
  static void uncontract(PartitionedHypergraph& partitioned_hg,
                         const Batch& batch,
                         gain_cache_t gain_cache) {
    switch(gain_cache.type) {
      case GainPolicy::cut:
        partitioned_hg.uncontract(batch, cast<CutGainCache>(gain_cache)); break;
      case GainPolicy::km1:
        partitioned_hg.uncontract(batch, cast<Km1GainCache>(gain_cache)); break;
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case GainPolicy::cut_for_graphs:
        partitioned_hg.uncontract(batch, cast<GraphCutGainCache>(gain_cache)); break;
      #endif
      case GainPolicy::none: break;
    }
  }

  template<typename PartitionedHypergraph, typename ParallelHyperedge>
  static void restoreSinglePinAndParallelNets(PartitionedHypergraph& partitioned_hg,
                                              const vec<ParallelHyperedge>& hes_to_restore,
                                              gain_cache_t gain_cache) {
    switch ( gain_cache.type ) {
      case GainPolicy::cut:
        partitioned_hg.restoreSinglePinAndParallelNets(hes_to_restore,
          cast<CutGainCache>(gain_cache)); break;
      case GainPolicy::km1:
        partitioned_hg.restoreSinglePinAndParallelNets(hes_to_restore,
          cast<Km1GainCache>(gain_cache)); break;
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case GainPolicy::cut_for_graphs:
        partitioned_hg.restoreSinglePinAndParallelNets(hes_to_restore,
          cast<GraphCutGainCache>(gain_cache)); break;
      #endif
      case GainPolicy::none: break;
    }
  }

  template<typename PartitionedHypergraph>
  static bool checkTrackedPartitionInformation(PartitionedHypergraph& partitioned_hg,
                                               gain_cache_t gain_cache) {
    switch ( gain_cache.type ) {
      case GainPolicy::cut:
        return partitioned_hg.checkTrackedPartitionInformation(
          cast<CutGainCache>(gain_cache));
      case GainPolicy::km1:
        return partitioned_hg.checkTrackedPartitionInformation(
          cast<Km1GainCache>(gain_cache));
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case GainPolicy::cut_for_graphs:
        return partitioned_hg.checkTrackedPartitionInformation(
          cast<GraphCutGainCache>(gain_cache));
      #endif
      case GainPolicy::none: return false;
    }
    return false;
  }

  template<typename GainCache>
  static GainCache& cast(gain_cache_t gain_cache) {
    if ( gain_cache.type != GainCache::TYPE ) {
      ERR("Cannot cast" << gain_cache.type << "to" << GainCache::TYPE);
    }
    return *reinterpret_cast<GainCache*>(gain_cache.gain_cache);
  }

 private:
  template<typename GainCache>
  static gain_cache_t constructGainCache() {
    return gain_cache_t { reinterpret_cast<gain_cache_s*>(new GainCache()), GainCache::TYPE };
  }
};

using GainCacheTypes = kahypar::meta::Typelist<Km1GainCache,
                                               CutGainCache
                                               ENABLE_GRAPHS(COMMA GraphCutGainCache)>;

#define INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_GAIN_CACHE(C)                                      \
  INSTANTIATE_CLASS_MACRO_WITH_TYPE_TRAITS_AND_OTHER_CLASS(C, Km1GainCache)                       \
  INSTANTIATE_CLASS_MACRO_WITH_TYPE_TRAITS_AND_OTHER_CLASS(C, CutGainCache)                       \
  ENABLE_GRAPHS(INSTANTIATE_CLASS_MACRO_WITH_TYPE_TRAITS_AND_OTHER_CLASS(C, GraphCutGainCache))

}  // namespace mt_kahypar

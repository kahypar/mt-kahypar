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
#include "mt-kahypar/partition/refinement/fm/gain_cache/km1_gain_cache.h"
#ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
#include "mt-kahypar/partition/refinement/fm/gain_cache/cut_gain_cache_for_graphs.h"
#endif
#include "mt-kahypar/partition/refinement/fm/gain_cache/do_nothing_gain_cache.h"
#include "mt-kahypar/macros.h"

namespace mt_kahypar {

struct gain_cache_s;
typedef struct  {
  gain_cache_s* gain_cache;
  FMGainCacheType type;
} gain_cache_t;

class AbstractGainCache {

 public:
  static gain_cache_t constructGainCache(const FMGainCacheType& type) {
    switch(type) {
      case FMGainCacheType::km1_gain_cache: return constructGainCache<Km1GainCache>();
      // TODO: replace this with cut gain cache once implementation is available
      case FMGainCacheType::cut_gain_cache: return constructGainCache<Km1GainCache>();
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case FMGainCacheType::cut_gain_cache_for_graphs: return constructGainCache<GraphCutGainCache>();
      #endif
      case FMGainCacheType::none: return constructGainCache<DoNothingGainCache>();
    }
    return constructGainCache<DoNothingGainCache>();
  }

  template<typename GainCache>
  static GainCache& cast(gain_cache_t gain_cache) {
    if ( gain_cache.type != GainCache::TYPE ) {
      ERR("Cannot cast" << gain_cache.type << "to" << GainCache::TYPE);
    }
    return *reinterpret_cast<GainCache*>(gain_cache.gain_cache);
  }

  static void deleteGainCache(gain_cache_t gain_cache) {
    if ( gain_cache.gain_cache ) {
      switch(gain_cache.type) {
        case FMGainCacheType::km1_gain_cache:
        case FMGainCacheType::cut_gain_cache:
          delete reinterpret_cast<Km1GainCache*>(gain_cache.gain_cache); break;
        #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
        case FMGainCacheType::cut_gain_cache_for_graphs:
          delete reinterpret_cast<GraphCutGainCache*>(gain_cache.gain_cache); break;
        #endif
        case FMGainCacheType::none:
          delete reinterpret_cast<DoNothingGainCache*>(gain_cache.gain_cache); break;
      }
    }
  }

 private:
  template<typename GainCache>
  static gain_cache_t constructGainCache() {
    return gain_cache_t { reinterpret_cast<gain_cache_s*>(new GainCache()), GainCache::TYPE };
  }
};

using GainCacheTypes = kahypar::meta::Typelist<Km1GainCache
                                               ENABLE_GRAPHS(COMMA GraphCutGainCache),
                                               DoNothingGainCache>;
using FMGainCacheTypes = kahypar::meta::Typelist<Km1GainCache
                                                 ENABLE_GRAPHS(COMMA GraphCutGainCache)>;

#define INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_GAIN_CACHE(C)                                      \
  INSTANTIATE_CLASS_MACRO_WITH_TYPE_TRAITS_AND_OTHER_CLASS(C, Km1GainCache)                       \
  ENABLE_GRAPHS(INSTANTIATE_CLASS_MACRO_WITH_TYPE_TRAITS_AND_OTHER_CLASS(C, GraphCutGainCache))   \
  INSTANTIATE_CLASS_MACRO_WITH_TYPE_TRAITS_AND_OTHER_CLASS(C, DoNothingGainCache)

#define INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_FM_GAIN_CACHE(C)                                   \
  INSTANTIATE_CLASS_MACRO_WITH_TYPE_TRAITS_AND_OTHER_CLASS(C, Km1GainCache)                       \
  ENABLE_GRAPHS(INSTANTIATE_CLASS_MACRO_WITH_TYPE_TRAITS_AND_OTHER_CLASS(C, GraphCutGainCache))

}  // namespace mt_kahypar

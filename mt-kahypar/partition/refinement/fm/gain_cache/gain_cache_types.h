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
#include "mt-kahypar/partition/refinement/fm/gain_cache/cut_gain_cache_for_graphs.h"
#include "mt-kahypar/partition/refinement/fm/gain_cache/do_nothing_gain_cache.h"
#include "mt-kahypar/macros.h"

namespace mt_kahypar {

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

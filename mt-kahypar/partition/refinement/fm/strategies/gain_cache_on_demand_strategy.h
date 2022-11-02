/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/partition/refinement/fm/fm_commons.h"
#include "gain_cache_strategy.h"

namespace mt_kahypar {

  class GainCacheOnDemandStrategy : public GainCacheStrategy {
  public:

    static constexpr bool maintain_gain_cache_between_rounds = false;

    GainCacheOnDemandStrategy(const Context& context,
                              HypernodeID numNodes,
                              FMSharedData& sharedData,
                              FMStats& runStats) :
            GainCacheStrategy(context, numNodes, sharedData, runStats),
            gainCacheInitMem(context.partition.k, 0)
    { }

    // conflicting signatures. derived does not have const qualifier for PHG. base has const. compiler doesn't complain, so probably fine.
    template<typename PHG>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void insertIntoPQ(PHG& phg, const HypernodeID v, const SearchID previous_search_of_v) {
      if (sharedData.nodeTracker.releasedMarker != previous_search_of_v) {
        // node is claimed for the first time in this fm round --> initialize gain cache entry
        phg.initializeGainCacheEntry(v, gainCacheInitMem);
      }
      GainCacheStrategy::insertIntoPQ(phg, v, previous_search_of_v);
    }

    void memoryConsumption(utils::MemoryTreeNode *parent) const {
      GainCacheStrategy::memoryConsumption(parent);
      parent->addChild("Initial Gain Comp", gainCacheInitMem.size() * sizeof(Gain));
    }

  private:
    vec<Gain> gainCacheInitMem;
  };


}
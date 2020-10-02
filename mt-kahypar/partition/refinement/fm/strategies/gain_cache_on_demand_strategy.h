/*******************************************************************************
 * This file is part of MT-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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
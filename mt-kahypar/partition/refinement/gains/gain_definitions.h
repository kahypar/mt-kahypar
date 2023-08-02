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

#include "kahypar-resources/meta/typelist.h"
#include "kahypar-resources/meta/policy_registry.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/partition/refinement/gains/km1/km1_gain_cache.h"
#include "mt-kahypar/partition/refinement/gains/km1/km1_rollback.h"
#include "mt-kahypar/partition/refinement/gains/km1/km1_gain_computation.h"
#include "mt-kahypar/partition/refinement/gains/km1/km1_attributed_gains.h"
#include "mt-kahypar/partition/refinement/gains/km1/km1_flow_network_construction.h"
#include "mt-kahypar/partition/refinement/gains/cut/cut_gain_cache.h"
#include "mt-kahypar/partition/refinement/gains/cut/cut_rollback.h"
#include "mt-kahypar/partition/refinement/gains/cut/cut_gain_computation.h"
#include "mt-kahypar/partition/refinement/gains/cut/cut_attributed_gains.h"
#include "mt-kahypar/partition/refinement/gains/cut/cut_flow_network_construction.h"
#ifdef KAHYPAR_ENABLE_SOED_METRIC
#include "mt-kahypar/partition/refinement/gains/soed/soed_attributed_gains.h"
#include "mt-kahypar/partition/refinement/gains/soed/soed_gain_computation.h"
#include "mt-kahypar/partition/refinement/gains/soed/soed_gain_cache.h"
#include "mt-kahypar/partition/refinement/gains/soed/soed_rollback.h"
#include "mt-kahypar/partition/refinement/gains/soed/soed_flow_network_construction.h"
#endif
#ifdef KAHYPAR_ENABLE_STEINER_TREE_METRIC
#include "mt-kahypar/partition/refinement/gains/steiner_tree/steiner_tree_attributed_gains.h"
#include "mt-kahypar/partition/refinement/gains/steiner_tree/steiner_tree_gain_computation.h"
#include "mt-kahypar/partition/refinement/gains/steiner_tree/steiner_tree_gain_cache.h"
#include "mt-kahypar/partition/refinement/gains/steiner_tree/steiner_tree_rollback.h"
#include "mt-kahypar/partition/refinement/gains/steiner_tree/steiner_tree_flow_network_construction.h"
#endif
#ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
#ifdef KAHYPAR_ENABLE_STEINER_TREE_METRIC
#include "mt-kahypar/partition/refinement/gains/steiner_tree_for_graphs/steiner_tree_attributed_gains_for_graphs.h"
#include "mt-kahypar/partition/refinement/gains/steiner_tree_for_graphs/steiner_tree_gain_computation_for_graphs.h"
#include "mt-kahypar/partition/refinement/gains/steiner_tree_for_graphs/steiner_tree_gain_cache_for_graphs.h"
#include "mt-kahypar/partition/refinement/gains/steiner_tree_for_graphs/steiner_tree_flow_network_construction_for_graphs.h"
#endif
#include "mt-kahypar/partition/refinement/gains/cut_for_graphs/cut_gain_cache_for_graphs.h"
#include "mt-kahypar/partition/refinement/gains/cut_for_graphs/cut_attributed_gains_for_graphs.h"
#endif
#include "mt-kahypar/macros.h"

namespace mt_kahypar {

struct Km1GainTypes : public kahypar::meta::PolicyBase {
  using GainComputation = Km1GainComputation;
  using AttributedGains = Km1AttributedGains;
  using GainCache = Km1GainCache;
  using DeltaGainCache = DeltaKm1GainCache;
  using Rollback = Km1Rollback;
  using FlowNetworkConstruction = Km1FlowNetworkConstruction;
};

struct CutGainTypes : public kahypar::meta::PolicyBase {
  using GainComputation = CutGainComputation;
  using AttributedGains = CutAttributedGains;
  using GainCache = CutGainCache;
  using DeltaGainCache = DeltaCutGainCache;
  using Rollback = CutRollback;
  using FlowNetworkConstruction = CutFlowNetworkConstruction;
};

#ifdef KAHYPAR_ENABLE_SOED_METRIC
struct SoedGainTypes : public kahypar::meta::PolicyBase {
  using GainComputation = SoedGainComputation;
  using AttributedGains = SoedAttributedGains;
  using GainCache = SoedGainCache;
  using DeltaGainCache = DeltaSoedGainCache;
  using Rollback = SoedRollback;
  using FlowNetworkConstruction = SoedFlowNetworkConstruction;
};
#endif

#ifdef KAHYPAR_ENABLE_STEINER_TREE_METRIC
struct SteinerTreeGainTypes : public kahypar::meta::PolicyBase {
  using GainComputation = SteinerTreeGainComputation;
  using AttributedGains = SteinerTreeAttributedGains;
  using GainCache = SteinerTreeGainCache;
  using DeltaGainCache = DeltaSteinerTreeGainCache;
  using Rollback = SteinerTreeRollback;
  using FlowNetworkConstruction = SteinerTreeFlowNetworkConstruction;
};
#endif

#ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
struct CutGainForGraphsTypes : public kahypar::meta::PolicyBase {
  using GainComputation = CutGainComputation;
  using AttributedGains = GraphCutAttributedGains;
  using GainCache = GraphCutGainCache;
  using DeltaGainCache = DeltaGraphCutGainCache;
  using Rollback = Km1Rollback;
  using FlowNetworkConstruction = CutFlowNetworkConstruction;
};

#ifdef KAHYPAR_ENABLE_STEINER_TREE_METRIC
struct SteinerTreeForGraphsTypes : public kahypar::meta::PolicyBase {
  using GainComputation = GraphSteinerTreeGainComputation;
  using AttributedGains = GraphSteinerTreeAttributedGains;
  using GainCache = GraphSteinerTreeGainCache;
  using DeltaGainCache = GraphDeltaSteinerTreeGainCache;
  using Rollback = SteinerTreeRollback;
  using FlowNetworkConstruction = GraphSteinerTreeFlowNetworkConstruction;
};
#endif
#endif


using GainTypes = kahypar::meta::Typelist<Km1GainTypes,
                                          CutGainTypes
                                          ENABLE_SOED(COMMA SoedGainTypes)
                                          ENABLE_STEINER_TREE(COMMA SteinerTreeGainTypes)
                                          ENABLE_GRAPHS(COMMA CutGainForGraphsTypes)
                                          ENABLE_GRAPHS(ENABLE_STEINER_TREE(COMMA SteinerTreeForGraphsTypes))>;

#define INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_GAIN_TYPES(C)                                                                 \
  INSTANTIATE_CLASS_MACRO_WITH_TYPE_TRAITS_AND_OTHER_CLASS(C, Km1GainTypes)                                                  \
  INSTANTIATE_CLASS_MACRO_WITH_TYPE_TRAITS_AND_OTHER_CLASS(C, CutGainTypes)                                                  \
  ENABLE_SOED(INSTANTIATE_CLASS_MACRO_WITH_TYPE_TRAITS_AND_OTHER_CLASS(C, SoedGainTypes))                                    \
  ENABLE_STEINER_TREE(INSTANTIATE_CLASS_MACRO_WITH_TYPE_TRAITS_AND_OTHER_CLASS(C, SteinerTreeGainTypes))                     \
  ENABLE_GRAPHS(INSTANTIATE_CLASS_MACRO_WITH_TYPE_TRAITS_AND_OTHER_CLASS(C, CutGainForGraphsTypes))                          \
  ENABLE_GRAPHS(ENABLE_STEINER_TREE(INSTANTIATE_CLASS_MACRO_WITH_TYPE_TRAITS_AND_OTHER_CLASS(C, SteinerTreeForGraphsTypes)))

}  // namespace mt_kahypar

/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@student.kit.edu>
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/coarsening/coarsening_commons.h"
#include "mt-kahypar/partition/coarsening/i_uncoarsener.h"
#include "mt-kahypar/partition/coarsening/uncoarsener_base.h"
#include "mt-kahypar/utils/progress_bar.h"

namespace mt_kahypar {

  class MultilevelUncoarsener : public IUncoarsener,
                                private UncoarsenerBase {

 public:
    MultilevelUncoarsener(Hypergraph& hypergraph,
                          const Context& context,
                          UncoarseningData& uncoarseningData) :
      UncoarsenerBase(hypergraph, context, uncoarseningData),
      _current_level(0),
      _num_levels(0),
      _block_ids(hypergraph.initialNumNodes(), PartIdType(kInvalidPartition)),
      _current_metrics(),
      _progress(hypergraph.initialNumNodes(), 0, false) { }

  MultilevelUncoarsener(const MultilevelUncoarsener&) = delete;
  MultilevelUncoarsener(MultilevelUncoarsener&&) = delete;
  MultilevelUncoarsener & operator= (const MultilevelUncoarsener &) = delete;
  MultilevelUncoarsener & operator= (MultilevelUncoarsener &&) = delete;

 private:
  void initializeImpl() override;

  bool isTopLevelImpl() const override;

  void projectToNextLevelAndRefineImpl() override;

  void rebalancingImpl() override;

  HypernodeID currentNumberOfNodesImpl() const override;

  PartitionedHypergraph&& movePartitionedHypergraphImpl() override;

  void refine(PartitionedHypergraph& partitioned_hypergraph, const double time_limit);

  int _current_level;
  int _num_levels;
  ds::Array<PartIdType> _block_ids;
  Metrics _current_metrics;
  utils::ProgressBar _progress;
};

}

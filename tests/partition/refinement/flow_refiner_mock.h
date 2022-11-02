/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
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


#include "kahypar/meta/registrar.h"

#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/partition/refinement/flows/i_flow_refiner.h"

namespace mt_kahypar {

using RefineFunc = std::function<MoveSequence(const PartitionedHypergraph&, const Subhypergraph&, const size_t)>;

class FlowRefinerMockControl {

  #define NOOP_REFINE_FUNC [] (const PartitionedHypergraph&, const Subhypergraph&, const size_t) { return MoveSequence { {}, 0 }; }

 public:
  FlowRefinerMockControl(const FlowRefinerMockControl&) = delete;
  FlowRefinerMockControl & operator= (const FlowRefinerMockControl &) = delete;

  FlowRefinerMockControl(FlowRefinerMockControl&&) = delete;
  FlowRefinerMockControl & operator= (FlowRefinerMockControl &&) = delete;

  static FlowRefinerMockControl& instance() {
    static FlowRefinerMockControl instance;
    return instance;
  }

 private:
  explicit FlowRefinerMockControl() :
    max_num_blocks(2),
    refine_func(NOOP_REFINE_FUNC) { }

 public:

  void reset() {
    max_num_blocks = 2;
    refine_func = NOOP_REFINE_FUNC;
  }

  PartitionID max_num_blocks;
  RefineFunc refine_func;
};

class FlowRefinerMock final : public IFlowRefiner {

 public:
  explicit FlowRefinerMock(const Hypergraph&,
                               const Context& context) :
    _context(context),
    _max_num_blocks(FlowRefinerMockControl::instance().max_num_blocks),
    _num_threads(0),
    _refine_func(FlowRefinerMockControl::instance().refine_func) { }

  FlowRefinerMock(const FlowRefinerMock&) = delete;
  FlowRefinerMock(FlowRefinerMock&&) = delete;
  FlowRefinerMock & operator= (const FlowRefinerMock &) = delete;
  FlowRefinerMock & operator= (FlowRefinerMock &&) = delete;

  virtual ~FlowRefinerMock() = default;

 protected:

 private:
  void initializeImpl(const PartitionedHypergraph&) { }

  MoveSequence refineImpl(const PartitionedHypergraph& phg,
                          const Subhypergraph& sub_hg,
                          const HighResClockTimepoint&) {
    return _refine_func(phg, sub_hg, _num_threads);
  }

  PartitionID maxNumberOfBlocksPerSearchImpl() const {
    return _max_num_blocks;
  }

  void setNumThreadsForSearchImpl(const size_t num_threads) {
    _num_threads = num_threads;
  }

  const Context& _context;
  const PartitionID _max_num_blocks;
  size_t _num_threads;
  RefineFunc _refine_func;
};

#define REGISTER_FLOW_REFINER(id, refiner)                                          \
  static kahypar::meta::Registrar<FlowRefinementFactory> register_ ## refiner(      \
    id,                                                                                 \
    [](const Hypergraph& hypergraph, const Context& context) -> IFlowRefiner* {     \
    return new refiner(hypergraph, context);                                            \
  })

REGISTER_FLOW_REFINER(FlowAlgorithm::mock, FlowRefinerMock);

}  // namespace mt_kahypar

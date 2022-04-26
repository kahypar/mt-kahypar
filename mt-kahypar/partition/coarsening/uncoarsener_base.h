/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@kit.edu>
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/coarsening/coarsening_commons.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/partition/metrics.h"


namespace mt_kahypar {

class UncoarsenerBase {

 protected:
  static constexpr bool debug = false;

 public:
  UncoarsenerBase(Hypergraph& hypergraph,
                          const Context& context,
                          UncoarseningData& uncoarseningData) :
          _hg(hypergraph),
          _context(context),
          _uncoarseningData(uncoarseningData) {}

  UncoarsenerBase(const UncoarsenerBase&) = delete;
  UncoarsenerBase(UncoarsenerBase&&) = delete;
  UncoarsenerBase & operator= (const UncoarsenerBase &) = delete;
  UncoarsenerBase & operator= (UncoarsenerBase &&) = delete;

  virtual ~UncoarsenerBase() = default;

  protected:
    Hypergraph& _hg;
    const Context& _context;
    UncoarseningData& _uncoarseningData;

  protected:

    double refinementTimeLimit(const Context& context, const double time) {
      if ( context.refinement.fm.time_limit_factor != std::numeric_limits<double>::max() ) {
        const double time_limit_factor = std::max(1.0,  context.refinement.fm.time_limit_factor * context.partition.k);
        return std::max(5.0, time_limit_factor * time);
      } else {
        return std::numeric_limits<double>::max();
      }
    }

  Metrics initialize(PartitionedHypergraph& phg) {
    Metrics m = { 0, 0, 0.0 };
    tbb::parallel_invoke([&] {
      m.cut = metrics::hyperedgeCut(phg);
    }, [&] {
      m.km1 = metrics::km1(phg);
    });
    m.imbalance = metrics::imbalance(phg, _context);

    int64_t num_nodes = phg.initialNumNodes();
    int64_t num_edges = Hypergraph::is_graph ? phg.initialNumEdges() / 2 : phg.initialNumEdges();
    utils::Stats::instance().add_stat("initial_num_nodes", num_nodes);
    utils::Stats::instance().add_stat("initial_num_edges", num_edges);
    utils::Stats::instance().add_stat("initial_cut", m.cut);
    utils::Stats::instance().add_stat("initial_km1", m.km1);
    utils::Stats::instance().add_stat("initial_imbalance", m.imbalance);
    return m;
  }

};
}

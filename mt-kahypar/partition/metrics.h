/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <kahypar/partition/metrics.h>
#include "mt-kahypar/partition/context.h"

namespace mt_kahypar::metrics {

    struct ThreadSafeMetrics {

    private:
        using Mode = kahypar::Mode;
        using Objective = kahypar::Objective;

        CAtomic<HyperedgeWeight> cut;
        CAtomic<HyperedgeWeight> km1;
        parallel::AtomicWrapper<double> imbalance;

    public:

        ThreadSafeMetrics() = default;
        ThreadSafeMetrics(HyperedgeWeight cut, HyperedgeWeight km1, double imbalance) : cut(cut), km1(km1), imbalance(imbalance) {}

        HyperedgeWeight loadCut() {
            return cut.load(std::memory_order_acquire);
        }

        HyperedgeWeight loadKm1() {
            return km1.load(std::memory_order_acquire);
        }

        double loadImbalance() {
            return imbalance.load(std::memory_order_acquire);
        }


        bool update_cut_strong(const HyperedgeWeight desired) {
            HyperedgeWeight expected = loadCut();
            return cut.compare_exchange_strong(expected,desired);
        }

        bool update_km1_strong(const HyperedgeWeight desired) {
            HyperedgeWeight expected = loadKm1();
            return km1.compare_exchange_strong(expected,desired);
        }

        bool update_imbalance_strong(const double desired) {
            double expected = loadImbalance();
            return imbalance.compare_exchange_strong(expected,desired);
        }

        void fetch_add(const HyperedgeWeight value, const Mode mode, const Objective objective) {
            if (mode == Mode::direct_kway) {
                switch (objective) {
                    case Objective::cut:
                        cut.fetch_add(value, std::memory_order_acq_rel);
                        break;
                    case Objective::km1:
                        km1.fetch_add(value, std::memory_order_acq_rel);
                        break;
                    default:
                        LOG << "Unknown Objective";
                        exit(-1);
                }
            } else if (mode == Mode::recursive_bisection) {
                // in recursive bisection, km1 is also optimized via the cut net metric
                cut.fetch_add(value, std::memory_order_acq_rel);
            }
        }

        HyperedgeWeight getMetric(const Mode mode, const Objective objective) {
            if (mode == Mode::direct_kway) {
                switch (objective) {
                    case Objective::cut: return cut.load(std::memory_order_acquire);
                    case Objective::km1: return km1.load(std::memory_order_acquire);
                    default:
                        LOG << "Unknown Objective";
                        exit(-1);
                }
            }
            ASSERT(mode == Mode::recursive_bisection);
            // in recursive bisection, km1 is also optimized via the cut net metric
            return cut.load(std::memory_order_acquire);
        }

        // ! Not thread-safe. To be used when a non-thread-safe version of this metrics object is required.
        kahypar::Metrics unsafeLoadMetrics() {
            return {cut.load(std::memory_order_acquire),
                    km1.load(std::memory_order_acquire),
                    imbalance.load(std::memory_order_acquire)};
        }

        // Not thread-safe. To be used when this thread-safe metrics object is supposed to be updated with values from a
        // non-threadsafe metrics object in a single threaded environment.
        void unsafeStoreMetrics(kahypar::Metrics& sequential) {
            cut.store(sequential.cut, std::memory_order_acq_rel);
            km1.store(sequential.km1, std::memory_order_acq_rel);
            imbalance.store(sequential.imbalance, std::memory_order_acq_rel);
        }
    };

HyperedgeWeight hyperedgeCut(const PartitionedHypergraph& hypergraph, bool parallel = true);

HyperedgeWeight km1(const PartitionedHypergraph& hypergraph, bool parallel = true);

HyperedgeWeight soed(const PartitionedHypergraph& hypergraph, bool parallel = true);

bool isBalanced(const PartitionedHypergraph& phg, const Context& context);

HyperedgeWeight objective(
        const PartitionedHypergraph& hg,
        const kahypar::Objective& objective,
        bool parallel = true);

double imbalance(const PartitionedHypergraph& hypergraph, const Context& context);


}  // namespace mt_kahypar

/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Nikolai Maas <nikolai.maas@student.kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <vector>
#include <fstream>

#include "tbb/parallel_reduce.h"

#include "mt-kahypar/definitions.h"


namespace mt_kahypar {
namespace utils {

template<typename T>
double parallel_stdev(const std::vector<T>& data, const double avg, const size_t n) {
    return std::sqrt(tbb::parallel_reduce(
            tbb::blocked_range<size_t>(0UL, data.size()), 0.0,
            [&](tbb::blocked_range<size_t>& range, double init) -> double {
            double tmp_stdev = init;
            for ( size_t i = range.begin(); i < range.end(); ++i ) {
                tmp_stdev += (data[i] - avg) * (data[i] - avg);
            }
            return tmp_stdev;
            }, std::plus<double>()) / ( n- 1 ));
}

template<typename T, typename F>
double parallel_weighted_stdev(const std::vector<T>& data, const double avg, const double n, F weight_fn) {
    return std::sqrt(tbb::parallel_reduce(
            tbb::blocked_range<size_t>(0UL, data.size()), 0.0,
            [&](tbb::blocked_range<size_t>& range, double init) -> double {
            double tmp_stdev = init;
            for ( size_t i = range.begin(); i < range.end(); ++i ) {
                tmp_stdev += weight_fn(i) * (data[i] - avg) * (data[i] - avg);
            }
            return tmp_stdev;
            }, std::plus<double>()) / ( n- 1 ));
}

template<typename T>
double parallel_avg(const std::vector<T>& data, const size_t n) {
    return tbb::parallel_reduce(
            tbb::blocked_range<size_t>(0UL, data.size()), 0.0,
            [&](tbb::blocked_range<size_t>& range, double init) -> double {
            double tmp_avg = init;
            for ( size_t i = range.begin(); i < range.end(); ++i ) {
                tmp_avg += static_cast<double>(data[i]);
            }
            return tmp_avg;
            }, std::plus<double>()) / static_cast<double>(n);
}

template<typename T, typename F>
double parallel_weighted_avg(const std::vector<T>& data, const double n, F weight_fn) {
    return tbb::parallel_reduce(
            tbb::blocked_range<size_t>(0UL, data.size()), 0.0,
            [&](tbb::blocked_range<size_t>& range, double init) -> double {
            double tmp_avg = init;
            for ( size_t i = range.begin(); i < range.end(); ++i ) {
                tmp_avg += weight_fn(i) * static_cast<double>(data[i]);
            }
            return tmp_avg;
            }, std::plus<double>()) / static_cast<double>(n);
}

static inline double avgHyperedgeDegree(const Hypergraph& hypergraph) {
    if (Hypergraph::is_graph) {
        return 2;
    }
    return static_cast<double>(hypergraph.initialNumPins()) / hypergraph.initialNumEdges();
}

static inline double avgHypernodeDegree(const Hypergraph& hypergraph) {
    return static_cast<double>(hypergraph.initialNumPins()) / hypergraph.initialNumNodes();
}

static void printDistributionStatsToCSV(const Hypergraph& hypergraph, const std::string& outfile) {
    ALWAYS_ASSERT(outfile != "");
    std::ofstream out(outfile.c_str());
    out << "degree, node-weight, incident-edge-weight, ratio" << std::endl;
    for (HypernodeID hn: hypergraph.nodes()) {
        const HyperedgeID degree = hypergraph.nodeDegree(hn);
        const HypernodeWeight node_weight = hypergraph.nodeWeight(hn);
        HyperedgeWeight edge_weight = 0;
        for (HyperedgeID he: hypergraph.incidentEdges(hn)) {
            edge_weight += hypergraph.edgeWeight(he);
        }
        const double ratio = static_cast<double>(edge_weight) / static_cast<double>(node_weight);
        out << degree << "," << node_weight << "," << edge_weight << "," << ratio << std::endl;
    }
}

#ifdef USE_GRAPH_PARTITIONER
static void outputGraphvizFile(const Hypergraph& hypergraph, const std::string& outfile) {
    ALWAYS_ASSERT(outfile != "");
    std::ofstream out(outfile.c_str());
    out.precision(3);

    out << "graph {" << std::endl;
    out << "fixedsize=true;" << std::endl << std::endl;

    for (HypernodeID hn: hypergraph.nodes()) {
        HypernodeWeight w = hypergraph.nodeWeight(hn);
        double root = std::sqrt(w) / 4; // std::round(100 * std::sqrt(w)) / 100.0;
        double penwidth = std::pow(w, 0.4);
        out << hn << " [label=\"\",width=" << root << ", penwidth=" << penwidth  << ",height=" << root << ",color=red];" << std::endl;
    }

    for (HyperedgeID e: hypergraph.edges()) {
        double w = std::sqrt(hypergraph.edgeWeight(e)) / 8;
        if (w > 0.01) {
            out << hypergraph.edgeSource(e) << "--" << hypergraph.edgeTarget(e)
                << " [weight=" << w << ", penwidth=" << w << "];" << std::endl;
        }
    }
    out << "}" << std::endl;
}
#endif

} // namespace utils
} // namespace mt_kahypar

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
#include "mt-kahypar/datastructures/separated_nodes.h"


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
double parallel_weighted_stdev(const std::vector<T>& data, const double avg, const double n, F weight_fn,
                               std::function<double(const T&)> mapping = [](auto x) { return x; }) {
    return std::sqrt(tbb::parallel_reduce(
            tbb::blocked_range<size_t>(0UL, data.size()), 0.0,
            [&](tbb::blocked_range<size_t>& range, double init) -> double {
            double tmp_stdev = init;
            for ( size_t i = range.begin(); i < range.end(); ++i ) {
                tmp_stdev += weight_fn(i) * (mapping(data[i]) - avg) * (mapping(data[i]) - avg);
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
double parallel_weighted_avg(const std::vector<T>& data, const double n, F weight_fn,
                             std::function<double(const T&)> mapping = [](auto x) { return x; }) {
    return tbb::parallel_reduce(
            tbb::blocked_range<size_t>(0UL, data.size()), 0.0,
            [&](tbb::blocked_range<size_t>& range, double init) -> double {
            double tmp_avg = init;
            for ( size_t i = range.begin(); i < range.end(); ++i ) {
                tmp_avg += weight_fn(i) * mapping(data[i]);
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

static inline bool isMeshGraph(const Hypergraph& graph) {
  const HypernodeID num_nodes = graph.initialNumNodes();
  const double avg_hn_degree = utils::avgHypernodeDegree(graph);
  std::vector<HyperedgeID> hn_degrees;
  hn_degrees.resize(graph.initialNumNodes());
  graph.doParallelForAllNodes([&](const HypernodeID& hn) {
    hn_degrees[hn] = graph.nodeDegree(hn);
  });
  const double stdev_hn_degree = utils::parallel_stdev(hn_degrees, avg_hn_degree, num_nodes);
  if (stdev_hn_degree > avg_hn_degree / 2) {
    return false;
  }

  // test whether 99.9th percentile hypernode degree is at most 4 times the average degree
  tbb::enumerable_thread_specific<size_t> num_high_degree_nodes(0);
  graph.doParallelForAllNodes([&](const HypernodeID& node) {
    if (graph.nodeDegree(node) > 4 * avg_hn_degree) {
      num_high_degree_nodes.local() += 1;
    }
  });
  return num_high_degree_nodes.combine(std::plus<>()) <= num_nodes / 1000;
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
namespace _private {
template<bool ColorBlocks>
struct Proxy { };

template<>
struct Proxy<false> {
  template<typename HypergraphT>
  static const char* getNodeColor(const HypergraphT&, HypernodeID, bool separated) {
    return separated ? "green": "blue";
  }

  template<typename HypergraphT>
  static bool isCut(const HypergraphT&, HyperedgeID) {
    return false;
  }

  template<typename HypergraphT>
  static PartitionID partID(const HypergraphT&, HypernodeID) {
    return kInvalidPartition;
  }

  template<typename HypergraphT>
  static PartitionID sepPartID(const HypergraphT&, HypernodeID) {
    return kInvalidPartition;
  }
};

template<>
struct Proxy<true> {
  template<typename HypergraphT>
  static const char* getNodeColor(const HypergraphT& hg, HypernodeID hn, bool) {
    // static const std::vector<const char*> colors{"red", "blue", "purple", "orange"};
    // ALWAYS_ASSERT(hg.partID(hn) != kInvalidPartition);
    // return colors.at(hg.partID(hn));
    ALWAYS_ASSERT(false);
    return "";
  }

  template<typename HypergraphT>
  static bool isCut(const HypergraphT& hg, HyperedgeID he) {
    return hg.partID(hg.edgeSource(he)) != hg.partID(hg.edgeTarget(he));
  }

  template<typename HypergraphT>
  static PartitionID partID(const HypergraphT& hg, HypernodeID hn) {
    ALWAYS_ASSERT(hg.partID(hn) != kInvalidPartition);
    return hg.partID(hn);
  }

  template<typename HypergraphT>
  static PartitionID sepPartID(const HypergraphT& hg, HypernodeID hn) {
    ALWAYS_ASSERT(hg.partID(hn) != kInvalidPartition);
    return hg.separatedPartID(hn);
  }
};

} // namespace _private

template<typename HypergraphT, bool ColorBlocks>
static const char* getNodeColor(const HypergraphT& hg, HypernodeID hn, bool separated = false) {
    return _private::Proxy<ColorBlocks>::getNodeColor(hg, hn, separated);
}

template<typename HypergraphT, bool ColorBlocks>
static bool isCut(const HypergraphT& hg, HyperedgeID he) {
    return _private::Proxy<ColorBlocks>::isCut(hg, he);
}

template<typename HypergraphT, bool ColorBlocks>
static PartitionID partID(const HypergraphT& hg, HypernodeID hn) {
    return _private::Proxy<ColorBlocks>::partID(hg, hn);
}

template<typename HypergraphT, bool ColorBlocks>
static PartitionID sepPartID(const HypergraphT& hg, HypernodeID hn) {
    return _private::Proxy<ColorBlocks>::sepPartID(hg, hn);
}

template<typename HypergraphT, bool ColorBlocks=HypergraphT::is_partitioned>
static void outputGraphvizFile(const HypergraphT& hypergraph, const std::string& outfile,
                               bool includeSeparated = false, const std::string& suffix = "") {
    // unused(hypergraph); unused(outfile); unused(includeSeparated); unused(suffix);
    ALWAYS_ASSERT(outfile != "");
    std::string file = outfile + suffix;
    std::ofstream out(file.c_str());
    out.precision(3);

    out << "graph {" << std::endl;
    out << "fixedsize=true;" << std::endl << std::endl;

    for (HypernodeID hn: hypergraph.nodes()) {
        HypernodeWeight w = hypergraph.nodeWeight(hn);
        double root = std::sqrt(w) / 4; // std::round(100 * std::sqrt(w)) / 100.0;
        double penwidth = std::pow(w, 0.4);
        const char* color = getNodeColor<HypergraphT, ColorBlocks>(hypergraph, hn);
        out << hn << " [label=\"\",width=" << root << ", penwidth=" << penwidth  << ",height="
            << root << ",color=" << color << "];" << std::endl;
    }
    if (includeSeparated && hypergraph.hasSeparatedNodes()) {
        const ds::SeparatedNodes& sn = hypergraph.separatedNodes().coarsest();
        for (HypernodeID sep: sn.nodes()) {
            HypernodeWeight w = sn.nodeWeight(sep);
            double root = std::sqrt(w) / 4;
            double penwidth = std::pow(w, 0.4);
            const char* color = getNodeColor<ds::SeparatedNodes, ColorBlocks>(sn, sep, true);
            out << sep + hypergraph.initialNumNodes() << " [label=\"\",width=" << root << ", penwidth=" << penwidth  << ",height="
                << root << ",color=" << color << "];" << std::endl;
        }
    }


    for (HyperedgeID e: hypergraph.edges()) {
        if (hypergraph.edgeTarget(e) >= hypergraph.edgeSource(e)) {
            // don't duplicate edges
            continue;
        }
        bool is_cut = isCut<HypergraphT, ColorBlocks>(hypergraph, e);
        double w = std::sqrt(hypergraph.edgeWeight(e)) / 8;
        if (w > 0.01) {
            out << hypergraph.edgeSource(e) << "--" << hypergraph.edgeTarget(e)
                << " [weight=" << w << ", penwidth=" << (is_cut ? 10 * w : w)
                << (is_cut ? ", color=green" : "")
                << "];" << std::endl;
        }
    }
    if (includeSeparated && hypergraph.hasSeparatedNodes()) {
        const ds::SeparatedNodes& sn = hypergraph.separatedNodes().coarsest();
        for (HypernodeID sep: sn.nodes()) {
            for (auto e: sn.inwardEdges(sep)) {
                bool is_cut = partID<HypergraphT, ColorBlocks>(hypergraph, e.target) != sepPartID<HypergraphT, ColorBlocks>(hypergraph, sep);
                double w = std::sqrt(e.weight) / 8;
                if (w > 0.01) {
                    out << sep + hypergraph.initialNumNodes() << "--" << e.target
                        << " [weight=" << w << ", penwidth=" << (is_cut ? 10 * w : w)
                        << (is_cut ? ", color=green" : "")
                        << "];" << std::endl;
                }
            }
        }
    }
    out << "}" << std::endl;
}
#endif

} // namespace utils
} // namespace mt_kahypar

#include <gtest/gtest.h>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/delete.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/datastructures/compressed_hypergraph_factory.h"

using ::testing::Test;
using namespace mt_kahypar;

namespace {

std::string write_temp_hgr(const std::string& name, const std::string& content) {
  std::string path = name; // relative to build dir when tests run
  std::ofstream out(path);
  out << content;
  out.close();
  return path;
}

std::vector<HypernodeID> pins_of(const ds::StaticHypergraph& hg, HyperedgeID e) {
  std::vector<HypernodeID> v;
  for (auto u : hg.pins(e)) v.push_back(u);
  return v;
}

std::vector<HypernodeID> pins_of(const ds::CompressedHypergraph& hg, HyperedgeID e) {
  std::vector<HypernodeID> v;
  for (auto u : hg.pins(e)) v.push_back(u);
  return v;
}

std::vector<HyperedgeID> incident_of(const ds::StaticHypergraph& hg, HypernodeID u) {
  std::vector<HyperedgeID> v;
  for (auto e : hg.incidentEdges(u)) v.push_back(e);
  return v;
}

std::vector<HyperedgeID> incident_of(const ds::CompressedHypergraph& hg, HypernodeID u) {
  std::vector<HyperedgeID> v;
  for (auto e : hg.incidentEdges(u)) v.push_back(e);
  return v;
}

} // namespace

// Helper to compare a static and compressed hypergraph thoroughly
static void assert_equality(const ds::StaticHypergraph& s, const ds::CompressedHypergraph& c) {
  // Global stats
  EXPECT_EQ(s.initialNumNodes(), c.initialNumNodes()) << "num nodes";
  EXPECT_EQ(s.initialNumEdges(), c.initialNumEdges()) << "num edges";
  EXPECT_EQ(s.totalWeight(),     c.totalWeight())     << "total weight";

  // Derived totals are more robust across loader implementations
  auto sum_edge_pins = [](const auto& hg) {
    size_t total = 0;
    for (auto e : hg.edges()) if (hg.edgeIsEnabled(e)) total += hg.edgeSize(e);
    return total;
  };
  auto sum_node_deg = [](const auto& hg) {
    size_t total = 0;
    for (auto u : hg.nodes()) if (hg.nodeIsEnabled(u)) total += hg.nodeDegree(u);
    return total;
  };
  const size_t s_pins = sum_edge_pins(s);
  const size_t c_pins = sum_edge_pins(c);
  const size_t s_deg  = sum_node_deg(s);
  const size_t c_deg  = sum_node_deg(c);
  if (s_pins != c_pins || s_deg != c_deg) {
    // Try to locate the first mismatching edge size or node degree
    const HyperedgeID m_chk = std::min<HyperedgeID>(s.initialNumEdges(), c.initialNumEdges());
    for (HyperedgeID e = 0; e < m_chk; ++e) {
      if (s.edgeIsEnabled(e) != c.edgeIsEnabled(e)) {
        ADD_FAILURE() << "edge enabled mismatch at e=" << e << ": s=" << s.edgeIsEnabled(e) << ", c=" << c.edgeIsEnabled(e);
        break;
      }
      if (!s.edgeIsEnabled(e)) continue;
      if (s.edgeSize(e) != c.edgeSize(e)) {
        ADD_FAILURE() << "edge size mismatch at e=" << e << ": s.size=" << s.edgeSize(e) << ", c.size=" << c.edgeSize(e);
        auto ps = pins_of(s, e);
        auto pc = pins_of(c, e);
        ADD_FAILURE() << "s pins (" << ps.size() << ") vs c pins (" << pc.size() << ")";
        break;
      }
    }
    const HypernodeID n_chk = std::min<HypernodeID>(s.initialNumNodes(), c.initialNumNodes());
    for (HypernodeID u = 0; u < n_chk; ++u) {
      if (s.nodeIsEnabled(u) != c.nodeIsEnabled(u)) {
        ADD_FAILURE() << "node enabled mismatch at u=" << u << ": s=" << s.nodeIsEnabled(u) << ", c=" << c.nodeIsEnabled(u);
        break;
      }
      if (!s.nodeIsEnabled(u)) continue;
      if (s.nodeDegree(u) != c.nodeDegree(u)) {
        ADD_FAILURE() << "node degree mismatch at u=" << u << ": s.deg=" << s.nodeDegree(u) << ", c.deg=" << c.nodeDegree(u);
        auto ia = incident_of(s, u);
        auto ib = incident_of(c, u);
        std::sort(ia.begin(), ia.end()); ia.erase(std::unique(ia.begin(), ia.end()), ia.end());
        std::sort(ib.begin(), ib.end()); ib.erase(std::unique(ib.begin(), ib.end()), ib.end());
        ADD_FAILURE() << "s incident size=" << ia.size() << ", c incident size=" << ib.size();
        break;
      }
    }
  }
  EXPECT_EQ(s_pins, c_pins) << "derived total pins differ";
  EXPECT_EQ(s_deg,  c_deg)  << "derived total degree differ";
  EXPECT_EQ(s_pins, s_deg)  << "static pins!=degree";
  EXPECT_EQ(c_pins, c_deg)  << "compressed pins!=degree";

  const HypernodeID n = s.initialNumNodes();
  const HyperedgeID m = s.initialNumEdges();

  // Nodes: enabled, weight, degree, incident sets
  for (HypernodeID u = 0; u < n; ++u) {
    EXPECT_EQ(s.nodeIsEnabled(u), c.nodeIsEnabled(u)) << "node enabled mismatch at u=" << u;
    if (!s.nodeIsEnabled(u)) continue;
    EXPECT_EQ(s.nodeWeight(u), c.nodeWeight(u)) << "node weight mismatch at u=" << u;
    EXPECT_EQ(s.nodeDegree(u), c.nodeDegree(u)) << "node degree mismatch at u=" << u;

    auto ia = incident_of(s, u);
    auto ib = incident_of(c, u);
    std::sort(ia.begin(), ia.end()); ia.erase(std::unique(ia.begin(), ia.end()), ia.end());
    std::sort(ib.begin(), ib.end()); ib.erase(std::unique(ib.begin(), ib.end()), ib.end());
    EXPECT_EQ(ia, ib) << "incident edges mismatch at u=" << u;
  }

  // Edges: enabled, weight, size, pins (in order)
  for (HyperedgeID e = 0; e < m; ++e) {
    EXPECT_EQ(s.edgeIsEnabled(e), c.edgeIsEnabled(e)) << "edge enabled mismatch at e=" << e;
    if (!s.edgeIsEnabled(e)) continue;
    EXPECT_EQ(s.edgeWeight(e), c.edgeWeight(e)) << "edge weight mismatch at e=" << e;
    EXPECT_EQ(s.edgeSize(e),   c.edgeSize(e))   << "edge size mismatch at e=" << e;

    auto ps = pins_of(s, e);
    auto pc = pins_of(c, e);
    ASSERT_EQ(ps.size(), pc.size()) << "pin count mismatch at e=" << e;
    for (size_t i = 0; i < ps.size(); ++i) {
      EXPECT_EQ(ps[i], pc[i]) << "pin mismatch at e=" << e << ", i=" << i;
    }
  }
}

// ---- Challenging graph generators (multi-byte varints for IDs and sizes) ----
namespace {
  // Build a large HGR content with:
  // - one very large edge with >127 pins (size header uses multi-byte varint)
  // - edges with large node IDs to force multi-byte varints for gaps
  // - a single-pin edge (to check loader removal)
  std::string make_large_hgr_content() {
    const int num_nodes = 5000;   // max node ID in file is <= num_nodes
    // We will define 5 edges
    // e0: large edge with 200 pins (IDs 200..399 inclusive -> 200 pins)
    // e1: unsorted with duplicates and large IDs
    // e2: single-pin edge
    // e3: small consecutive block to keep variety
    // e4: three large IDs far apart
    std::string content;
    content.reserve(4096);
    content += "5 ";
    content += std::to_string(num_nodes);
    content += "\n";

    // e0: 200 pins -> 200..399 (1-based already fits format used in previous tests)
    for (int i = 200; i < 400; ++i) {
      content += std::to_string(i);
      if (i + 1 < 400) content += ' ';
    }
    content += "\n";

    // e1: duplicates + unsorted including large IDs
    content += "300 172 172 1 4096 3500 4096\n";

    // e2: single pin
    content += "5\n";

    // e3: small block
    content += "200 201 202 203 204\n";

    // e4: widely spaced large IDs (first gap > 127)
    content += "1000 1300 2600\n";

    return content;
  }

  // Build in-memory edges for factory/contract tests mirroring the above shape.
  void build_large_edges(parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>& edges,
                         HypernodeID& num_nodes_out) {
    const HypernodeID num_nodes = 5000;
    const HyperedgeID num_edges = 5;
    edges.resize(num_edges);

    // e0: 200 pins [200..399] (edge ID 0-based IDs for internal APIs)
    edges[0].reserve(200);
    for (HypernodeID u = 200; u < 400; ++u) edges[0].push_back(u);

    // e1: duplicates + unsorted including large IDs
    edges[1] = { 300, 172, 172, 1, 4096, 3500, 4096 };

    // e2: single pin
    edges[2] = { 5 };

    // e3: small block
    edges[3] = { 200, 201, 202, 203, 204 };

    // e4: widely spaced large IDs
    edges[4] = { 1000, 1300, 2600 };

    num_nodes_out = num_nodes;
  }
}

TEST(CompressedHypergraphIO, MatchesStaticHypergraphIO) {
  const std::string content = make_large_hgr_content();

  const std::string path = write_temp_hgr("tmp_comp_io_test.hgr", content);

  // Read as static
  mt_kahypar_hypergraph_t h1 = io::readInputFile(path, PresetType::default_preset,
                                                 InstanceType::hypergraph,
                                                 FileFormat::hMetis, false /*stable*/);
  const auto& s = utils::cast<ds::StaticHypergraph>(h1);

  // Read as compressed
  mt_kahypar_hypergraph_t h2 = io::readInputFile(path, PresetType::default_preset,
                                                 InstanceType::compressed_hypergraph,
                                                 FileFormat::hMetis, false /*stable*/);
  const auto& c = utils::cast<ds::CompressedHypergraph>(h2);

  assert_equality(s, c);

  // Cleanup
  utils::delete_hypergraph(h1);
  utils::delete_hypergraph(h2);
}

TEST(CompressedHypergraphRemoval, RemoveEdgeAndRemoveLargeEdgeMatch) {
  const std::string content = make_large_hgr_content();

  const std::string path = write_temp_hgr("tmp_remove_edge_test.hgr", content);

  // Load both representations
  mt_kahypar_hypergraph_t h1 = io::readInputFile(path, PresetType::default_preset,
                                                 InstanceType::hypergraph,
                                                 FileFormat::hMetis, false /*stable*/);
  mt_kahypar_hypergraph_t h2 = io::readInputFile(path, PresetType::default_preset,
                                                 InstanceType::compressed_hypergraph,
                                                 FileFormat::hMetis, false /*stable*/);
  auto& s = const_cast<ds::StaticHypergraph&>(utils::cast<ds::StaticHypergraph>(h1));
  auto& c = const_cast<ds::CompressedHypergraph&>(utils::cast<ds::CompressedHypergraph>(h2));

  // After loader removes the single-pin edge, we should have 4 remaining
  ASSERT_EQ(s.initialNumEdges(), 4);
  ASSERT_EQ(c.initialNumEdges(), 4);

  // Sanity: edges should be enabled in both
  ASSERT_TRUE(s.edgeIsEnabled(0));
  ASSERT_TRUE(s.edgeIsEnabled(1));
  ASSERT_TRUE(s.edgeIsEnabled(2));
  ASSERT_TRUE(c.edgeIsEnabled(0));
  ASSERT_TRUE(c.edgeIsEnabled(1));
  ASSERT_TRUE(c.edgeIsEnabled(2));

  // Remove edge 0 using removeEdge in both (this is the large edge)
  s.removeEdge(0);
  c.removeEdge(0);
  EXPECT_FALSE(s.edgeIsEnabled(0));
  EXPECT_FALSE(c.edgeIsEnabled(0));

  // Remove edge 1 using removeLargeEdge in both
  s.removeLargeEdge(1);
  c.removeLargeEdge(1);
  EXPECT_FALSE(s.edgeIsEnabled(1));
  EXPECT_FALSE(c.edgeIsEnabled(1));

  // Full-graph equality after removals
  assert_equality(s, c);

  utils::delete_hypergraph(h1);
  utils::delete_hypergraph(h2);
}

TEST(CompressedHypergraphRestoration, RestoreEdgeAndRestoreLargeEdgeWork) {
  const std::string content = make_large_hgr_content();

  const std::string path = write_temp_hgr("tmp_restore_edge_test.hgr", content);

  // Load both representations
  mt_kahypar_hypergraph_t h1 = io::readInputFile(path, PresetType::default_preset,
                                                 InstanceType::hypergraph,
                                                 FileFormat::hMetis, false /*stable*/);
  mt_kahypar_hypergraph_t h2 = io::readInputFile(path, PresetType::default_preset,
                                                 InstanceType::compressed_hypergraph,
                                                 FileFormat::hMetis, false /*stable*/);
  auto& s = const_cast<ds::StaticHypergraph&>(utils::cast<ds::StaticHypergraph>(h1));
  auto& c = const_cast<ds::CompressedHypergraph&>(utils::cast<ds::CompressedHypergraph>(h2));

  assert_equality(s, c);

  // Remove two edges via different paths
  s.removeEdge(0);
  c.removeEdge(0);
  s.removeLargeEdge(1);
  c.removeLargeEdge(1);

  assert_equality(s, c);

  // Restore both back
  s.restoreLargeEdge(0);
  c.restoreLargeEdge(0);
  s.restoreLargeEdge(1);
  c.restoreLargeEdge(1);

  // After restoration, graphs must match again
  assert_equality(s, c);

  utils::delete_hypergraph(h1);
  utils::delete_hypergraph(h2);
}

TEST(CompressedHypergraphContract, ContractsLikeStatic) {
  // Build challenging in-memory graph mirroring HGR
  parallel::scalable_vector<parallel::scalable_vector<HypernodeID>> edges;
  HypernodeID num_nodes;
  build_large_edges(edges, num_nodes);

  const HyperedgeID num_edges = edges.size();
  // Normalize pins per edge to match loader semantics (sorted unique)
  for (auto& e : edges) { std::sort(e.begin(), e.end()); e.erase(std::unique(e.begin(), e.end()), e.end()); }
  std::vector<HyperedgeWeight> he_w(num_edges, 1);
  std::vector<HypernodeWeight> hn_w(num_nodes, 1);

  ds::StaticHypergraph s = ds::StaticHypergraphFactory::construct(num_nodes, num_edges, edges,
                                                                  he_w.data(), hn_w.data(), false);
  ds::CompressedHypergraph c = ds::CompressedHypergraphFactory::construct(num_nodes, num_edges, edges,
                                                                          he_w.data(), hn_w.data(), false);

  // Partition nodes into 8 blocks by ranges to induce non-trivial contraction
  parallel::scalable_vector<HypernodeID> communities(num_nodes);
  for (HypernodeID u = 0; u < num_nodes; ++u) {
    communities[u] = static_cast<HypernodeID>(u / 700); // ~8 groups up to 5000
  }

  auto s_coarse = s.contract(communities, /*deterministic=*/false);
  auto c_coarse = c.contract(communities, /*deterministic=*/false);

  assert_equality(s_coarse, c_coarse);
}

TEST(CompressedHypergraphFactoryConstruct, MatchesStaticFactoryConstruct) {
  parallel::scalable_vector<parallel::scalable_vector<HypernodeID>> edges;
  HypernodeID num_nodes;
  build_large_edges(edges, num_nodes);
  const HyperedgeID num_edges = edges.size();

  // Normalize pins per edge to match loader semantics (sorted unique)
  for (auto& e : edges) { std::sort(e.begin(), e.end()); e.erase(std::unique(e.begin(), e.end()), e.end()); }

  std::vector<HyperedgeWeight> he_w(num_edges, 1);
  std::vector<HypernodeWeight> hn_w(num_nodes, 1);

  ds::StaticHypergraph s = ds::StaticHypergraphFactory::construct(num_nodes, num_edges, edges,
                                                                  he_w.data(), hn_w.data(), false);
  ds::CompressedHypergraph c = ds::CompressedHypergraphFactory::construct(num_nodes, num_edges, edges,
                                                                          he_w.data(), hn_w.data(), false);

  assert_equality(s, c);
}

// --- Two-level weights specific tests ---

TEST(CompressedHypergraphTwoLevelWeights, DefaultIsOneAndLazy) {
  // Build edges and construct compressed hypergraph WITHOUT explicit weights.
  parallel::scalable_vector<parallel::scalable_vector<HypernodeID>> edges;
  HypernodeID num_nodes;
  build_large_edges(edges, num_nodes);
  const HyperedgeID num_edges = edges.size();

  // Normalize pins per edge to match loader semantics (sorted unique)
  for (auto& e : edges) { std::sort(e.begin(), e.end()); e.erase(std::unique(e.begin(), e.end()), e.end()); }

  ds::CompressedHypergraph c = ds::CompressedHypergraphFactory::construct(
      num_nodes, num_edges, edges, /*edge_weights=*/nullptr, /*node_weights=*/nullptr, /*stable=*/false);

  // All node and edge weights should read as 1; total weight equals number of nodes
  EXPECT_EQ(static_cast<HypernodeWeight>(num_nodes), c.totalWeight()) << "total node weight should be V when no weights provided";
  for (HypernodeID u = 0; u < num_nodes; ++u) {
    EXPECT_EQ(1, c.nodeWeight(u)) << "default node weight must be 1 at u=" << u;
  }
  for (HyperedgeID e = 0; e < num_edges; ++e) {
    EXPECT_EQ(1, c.edgeWeight(e)) << "default edge weight must be 1 at e=" << e;
  }
}

TEST(CompressedHypergraphTwoLevelWeights, SupportsLargeAndSparseValuesAndContract) {
  // Build edges and assign mostly-1 weights with a few large values to exercise overflow path.
  parallel::scalable_vector<parallel::scalable_vector<HypernodeID>> edges;
  HypernodeID num_nodes;
  build_large_edges(edges, num_nodes);
  const HyperedgeID num_edges = edges.size();

  // Normalize pins per edge to match loader semantics (sorted unique)
  for (auto& e : edges) { std::sort(e.begin(), e.end()); e.erase(std::unique(e.begin(), e.end()), e.end()); }

  std::vector<HyperedgeWeight> he_w(num_edges, 1);
  std::vector<HypernodeWeight> hn_w(num_nodes, 1);
  // Large/sparse values
  he_w[1] = 300;        // > 255 to exceed uint8_t base
  he_w[3] = 10000;      // much larger
  hn_w[0] = 500;        // > 255
  hn_w[2600] = 70000;   // large value

  ds::StaticHypergraph s = ds::StaticHypergraphFactory::construct(num_nodes, num_edges, edges,
                                                                  he_w.data(), hn_w.data(), false);
  ds::CompressedHypergraph c = ds::CompressedHypergraphFactory::construct(num_nodes, num_edges, edges,
                                                                          he_w.data(), hn_w.data(), false);

  // The two representations must agree on weights and structure
  assert_equality(s, c);

  // Also validate that contraction preserves weight semantics identically in both reps
  parallel::scalable_vector<HypernodeID> communities(num_nodes);
  for (HypernodeID u = 0; u < num_nodes; ++u) communities[u] = static_cast<HypernodeID>(u / 700);
  auto s_coarse = s.contract(communities, /*deterministic=*/false);
  auto c_coarse = c.contract(communities, /*deterministic=*/false);
  assert_equality(s_coarse, c_coarse);
}

TEST(CompressedHypergraphTwoLevelWeights, SettersUpdateValuesCorrectly) {
  // Start with no explicit weights, then set a few nodes/edges to large values
  parallel::scalable_vector<parallel::scalable_vector<HypernodeID>> edges;
  HypernodeID num_nodes;
  build_large_edges(edges, num_nodes);
  const HyperedgeID num_edges = edges.size();
  for (auto& e : edges) { std::sort(e.begin(), e.end()); e.erase(std::unique(e.begin(), e.end()), e.end()); }

  ds::CompressedHypergraph c = ds::CompressedHypergraphFactory::construct(
      num_nodes, num_edges, edges, /*edge_weights=*/nullptr, /*node_weights=*/nullptr, /*stable=*/false);

  // Verify defaults
  EXPECT_EQ(static_cast<HypernodeWeight>(num_nodes), c.totalWeight());
  EXPECT_EQ(1, c.nodeWeight(0));
  EXPECT_EQ(1, c.edgeWeight(0));

  // Update a few weights, including very large values
  c.setNodeWeight(0, static_cast<HypernodeWeight>(500));
  c.setNodeWeight(2600, static_cast<HypernodeWeight>(70000));
  c.setEdgeWeight(1, static_cast<HyperedgeWeight>(300));
  c.setEdgeWeight(3, static_cast<HyperedgeWeight>(10000));

  // Read back
  EXPECT_EQ(500, c.nodeWeight(0));
  EXPECT_EQ(70000, c.nodeWeight(2600));
  EXPECT_EQ(300, c.edgeWeight(1));
  EXPECT_EQ(10000, c.edgeWeight(3));

  // Total node weight should reflect updates: started as num_nodes (all ones), then +499 and +69999
  const HypernodeWeight expected_total = static_cast<HypernodeWeight>(num_nodes)
                                       + static_cast<HypernodeWeight>(499)
                                       + static_cast<HypernodeWeight>(69999);
  EXPECT_EQ(expected_total, c.totalWeight());
}

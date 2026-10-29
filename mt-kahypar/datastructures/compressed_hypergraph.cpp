/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "compressed_hypergraph.h"

#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/datastructures/concurrent_bucket_map.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/memory_tree.h"

#include <tbb/parallel_reduce.h>
#include <tbb/parallel_sort.h>
#include <tbb/enumerable_thread_specific.h>

#include <unordered_map>
#include <algorithm>
#include <numeric>

namespace mt_kahypar::ds {

  // Raw varint helpers for fixed-position encoding to a byte buffer
  static inline size_t varint_size_u64(uint64_t v) {
    size_t n = 1;
    while (v >= 0x80) { v >>= 7; ++n; }
    return n;
  }

  static inline size_t encode_varint_to_ptr(uint64_t v, uint8_t* out) {
    size_t i = 0;
    while (v >= 0x80) {
      out[i++] = static_cast<uint8_t>((v & 0x7F) | 0x80);
      v >>= 7;
    }
    out[i++] = static_cast<uint8_t>(v & 0x7F);
    return i;
  }


  /*!
  * This struct is used during multilevel coarsening to efficiently
  * detect parallel hyperedges.
  */
  struct ContractedHyperedgeInformation {
    HyperedgeID he = kInvalidHyperedge;
    size_t hash = kEdgeHashSeed;
    size_t size = std::numeric_limits<size_t>::max();
    bool valid = false;
  };

  /*!
   * Contracts a given community structure. All vertices with the same label
   * are collapsed into the same vertex. The resulting single-pin and parallel
   * hyperedges are removed from the contracted graph. The function returns
   * the contracted hypergraph and a mapping which specifies a mapping from
   * community label (given in 'communities') to a vertex in the coarse hypergraph.
   *
   * \param communities Community structure that should be contracted
   */
  CompressedHypergraph CompressedHypergraph::contract(parallel::scalable_vector<HypernodeID>& communities, bool /* deterministic */) {
    // Stage 0: Preconditions and helpers
    ASSERT(communities.size() == _num_hypernodes);

    auto map_to_coarse = [&](const HypernodeID hn) -> HypernodeID {
      ASSERT(hn < communities.size());
      return communities[hn];
    };

    // #################### STAGE 1 ####################
    // Densify community IDs and compute number of coarse nodes
    // We assume communities[hn] in [0, _num_hypernodes), but only enabled nodes count.
    std::vector<uint8_t> present(_num_hypernodes, 0);
    for (HypernodeID hn = 0; hn < _num_hypernodes; ++hn) {
      if (nodeIsEnabled(hn)) {
        const HypernodeID cid = communities[hn];
        ASSERT(cid < _num_hypernodes);
        present[cid] = 1;
      } else {
        communities[hn] = kInvalidHypernode;
      }
    }

  std::vector<HypernodeID> label_to_coarse(_num_hypernodes, kInvalidHypernode);
    HypernodeID next_coarse = 0;
    for (HypernodeID lab = 0; lab < _num_hypernodes; ++lab) {
      if (present[lab]) label_to_coarse[lab] = next_coarse++;
    }
    const HypernodeID num_coarse_nodes = next_coarse;

    // Remap communities to dense IDs
    for (HypernodeID hn = 0; hn < _num_hypernodes; ++hn) {
      if (communities[hn] != kInvalidHypernode) {
        communities[hn] = label_to_coarse[communities[hn]];
      }
    }

    // Release temporary label maps early
    present.clear();
    present.shrink_to_fit();
    label_to_coarse.clear();
    label_to_coarse.shrink_to_fit();

    // Aggregate coarse node weights
    std::vector<HypernodeWeight> coarse_node_weight(num_coarse_nodes, 0);
    for (HypernodeID hn = 0; hn < _num_hypernodes; ++hn) {
      const HypernodeID ch = communities[hn];
      if (ch != kInvalidHypernode) {
        coarse_node_weight[ch] += nodeWeight(hn);
      }
    }

    // #################### STAGE 2 ####################
    // Deduplicate coarse edges without storing per-edge byte slices.
    // We collect candidate edges as IDs with lightweight metadata, sort by (hash, size),
    // and tie-break with on-the-fly lexicographic comparison of coarse pins.

    auto build_coarse_pins = [&](HyperedgeID he, std::vector<HypernodeID>& out) {
      auto HE = hyperedge(he);
      const size_t begin = HE.firstEntry();
      const size_t end   = HE.firstInvalidEntry();
      size_t pos = begin;
      HypernodeID prev = 0;
      const size_t esize = static_cast<size_t>(HE.size());
      out.clear();
      out.reserve(esize);
      size_t emitted = 0;
      while (pos < end && emitted < esize) {
        const uint64_t gap = decode_varint_bounded(_compressed_incidence_array, pos, end);
        if (emitted > 0 && gap == 0) { continue; }
        const HypernodeID v = static_cast<HypernodeID>(prev + gap);
        prev = v;
        ++emitted;
        const HypernodeID ch = map_to_coarse(v);
        if (ch != kInvalidHypernode) out.push_back(ch);
      }
      if (!out.empty()) {
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
      }
    };

    auto hash_pins = [&](const std::vector<HypernodeID>& v) -> uint64_t {
      // 64-bit FNV-1a style mix
      uint64_t h = 1469598103934665603ull;
      for (HypernodeID x : v) {
        h ^= static_cast<uint64_t>(x) + 0x9e3779b97f4a7c15ULL;
        h *= 1099511628211ull;
      }
      return h;
    };

    // Build candidate list and metadata in parallel with thread-local buffers
    struct LocalMeta {
      std::vector<HyperedgeID> cand;
      std::vector<uint64_t> chash;
      std::vector<uint32_t> csize;
      std::vector<uint32_t> cfp;
    };
    tbb::enumerable_thread_specific<LocalMeta> locals;
    tbb::enumerable_thread_specific<std::vector<HypernodeID>> tls_pins;
    tbb::parallel_for(HyperedgeID(0), _num_hyperedges, [&](const HyperedgeID he) {
      if (!edgeIsEnabled(he)) return;
      auto& pins = tls_pins.local();
      pins.clear();
      pins.reserve(32);
      build_coarse_pins(he, pins);
      if (pins.size() <= 1) return;
      auto& lm = locals.local();
      lm.cand.push_back(he);
      lm.chash.push_back(hash_pins(pins));
      lm.csize.push_back(static_cast<uint32_t>(pins.size()));
      // 32-bit FNV-1a fingerprint as a cheap lex locality key
      uint32_t fp = 2166136261u;
      for (HypernodeID v : pins) { fp ^= static_cast<uint32_t>(v); fp *= 16777619u; }
      lm.cfp.push_back(fp);
    });
    // Merge locals deterministically by thread enumeration order
    std::vector<HyperedgeID> cand;
    std::vector<uint64_t> chash;
    std::vector<uint32_t> csize;
    std::vector<uint32_t> cfp;
    size_t total = 0;
    locals.combine_each([&](const LocalMeta& lm) { total += lm.cand.size(); });
    cand.reserve(total); chash.reserve(total); csize.reserve(total); cfp.reserve(total);
    locals.combine_each([&](LocalMeta& lm) {
      cand.insert(cand.end(), lm.cand.begin(), lm.cand.end());
      chash.insert(chash.end(), lm.chash.begin(), lm.chash.end());
      csize.insert(csize.end(), lm.csize.begin(), lm.csize.end());
      cfp.insert(cfp.end(), lm.cfp.begin(), lm.cfp.end());
      // free local buffers
      lm.cand.clear(); lm.chash.clear(); lm.csize.clear(); lm.cfp.clear();
      lm.cand.shrink_to_fit(); lm.chash.shrink_to_fit(); lm.csize.shrink_to_fit(); lm.cfp.shrink_to_fit();
    });

    std::vector<size_t> perm(cand.size());
    std::iota(perm.begin(), perm.end(), 0);
    // Comparator avoids decoding; equality is resolved in the dedup pass
    auto cmp_coarse = [&](size_t ia, size_t ib) {
      const uint64_t ha = chash[ia], hb = chash[ib];
      if (ha != hb) return ha < hb;
      if (csize[ia] != csize[ib]) return csize[ia] < csize[ib];
      if (cfp[ia] != cfp[ib]) return cfp[ia] < cfp[ib];
      return cand[ia] < cand[ib];
    };
    tbb::parallel_sort(perm.begin(), perm.end(), cmp_coarse);

    // Build deduplicated representative list, weights, and encoded lengths
    std::vector<HyperedgeID> uniq_rep; uniq_rep.reserve(perm.size());
    std::vector<HyperedgeWeight> uniq_weight; uniq_weight.reserve(perm.size());
    std::vector<uint32_t> uniq_len; uniq_len.reserve(perm.size());
    std::vector<uint32_t> uniq_pins; uniq_pins.reserve(perm.size());

    size_t i = 0; std::vector<HypernodeID> a, b, cached;
    bool cached_valid = false; size_t cached_perm_index = 0;
    while (i < perm.size()) {
      const size_t idx = perm[i];
      const HyperedgeID rep = cand[idx];
      // Reuse cached pins if we already decoded them for the next representative
      if (cached_valid && cached_perm_index == idx) { a.swap(cached); cached_valid = false; }
      else { build_coarse_pins(rep, a); }
      const uint32_t pins_cnt = static_cast<uint32_t>(a.size());
      // Compute encoded length for representative
      size_t elen = varint_size_u64(pins_cnt);
      HypernodeID pprev = 0;
      for (HypernodeID v : a) { elen += varint_size_u64(static_cast<uint64_t>(v - pprev)); pprev = v; }
      HyperedgeWeight wsum = edgeWeight(rep);
      size_t j = i + 1;
      while (j < perm.size()) {
        const size_t idx2 = perm[j];
        if (chash[idx2] != chash[idx]) break;
        if (csize[idx2] != pins_cnt) break;
        // Fast pre-filter: skip decoding on fingerprint mismatch
        if (cfp[idx2] != cfp[idx]) { break; }
        build_coarse_pins(cand[idx2], b);
        if (b != a) {
          // carry b to next representative to avoid re-decoding
          cached.swap(b);
          cached_perm_index = idx2;
          cached_valid = true;
          break;
        }
        wsum += edgeWeight(cand[idx2]);
        ++j;
      }
      uniq_rep.push_back(rep);
      uniq_weight.push_back(wsum);
      uniq_len.push_back(static_cast<uint32_t>(elen));
      uniq_pins.push_back(pins_cnt);
      i = j;
    }

    // We can release candidate metadata now
    cand.clear(); cand.shrink_to_fit();
    chash.clear(); chash.shrink_to_fit();
    csize.clear(); csize.shrink_to_fit();
    perm.clear(); perm.shrink_to_fit();

    // Pre-compute total bytes for incidence encoding
    size_t total_inc_bytes = 0;
    for (uint32_t l : uniq_len) total_inc_bytes += l;

    // #################### STAGE 3 ####################
    // Build coarse hypergraph with compressed storage
    CompressedHypergraph hypergraph;
    hypergraph._num_hypernodes = num_coarse_nodes;
  hypergraph._num_hyperedges = static_cast<HyperedgeID>(uniq_rep.size());
    hypergraph._num_removed_hyperedges = 0;
    hypergraph._max_edge_size = 0;
    hypergraph._num_pins = 0;
    hypergraph._total_degree = 0;
    hypergraph._total_weight = 0;

    // Init CSR containers
    hypergraph._hypernode_offsets.assign(num_coarse_nodes, 0);
    hypergraph._hypernode_enabled.assign(num_coarse_nodes, 1);
    hypergraph._hyperedge_offsets.assign(hypergraph._num_hyperedges, 0);
    hypergraph._hyperedge_enabled.assign(hypergraph._num_hyperedges, 1);
    hypergraph._compressed_incidence_array.clear();
    if (total_inc_bytes > 0) hypergraph._compressed_incidence_array.reserve(total_inc_bytes);
    hypergraph._compressed_incident_nets.clear();

    // Set node weights
    hypergraph._hypernode_weights.ensure_initialized(num_coarse_nodes);
    for (HypernodeID u = 0; u < num_coarse_nodes; ++u) {
      hypergraph._hypernode_weights[u] = coarse_node_weight[u];
      // firstEntry/compressed sizes set later
    }
    // Coarse node weights no longer needed
    coarse_node_weight.clear();
    coarse_node_weight.shrink_to_fit();

    // Encode coarse edges (pins) directly into final buffer
    // Pre-size incidence array and compute offsets by prefix sum of uniq_len
    size_t pins_bytes_offset = 0;
    if (total_inc_bytes > 0) hypergraph._compressed_incidence_array.resize(total_inc_bytes);
    for (HyperedgeID he = 0; he < hypergraph._num_hyperedges; ++he) {
      hypergraph._hyperedge_offsets[he] = pins_bytes_offset;
      pins_bytes_offset += uniq_len[he];
    }
    if (hypergraph._hyperedge_weights.empty()) hypergraph._hyperedge_weights.ensure_initialized(hypergraph._num_hyperedges);
    // Assign weights sequentially to avoid potential contention in TwoLevelVector
    for (HyperedgeID he = 0; he < hypergraph._num_hyperedges; ++he) {
      hypergraph._hyperedge_weights[he] = uniq_weight[he];
    }
    // Encode pins in parallel into pre-sized buffer
    uint8_t* base_ptr = hypergraph._compressed_incidence_array.data();
    tbb::parallel_for(HyperedgeID(0), hypergraph._num_hyperedges, [&](const HyperedgeID he) {
      std::vector<HypernodeID> pins;
      pins.reserve(32);
      build_coarse_pins(uniq_rep[he], pins);
      uint8_t* ptr = base_ptr + hypergraph._hyperedge_offsets[he];
      ptr += encode_varint_to_ptr(static_cast<uint64_t>(pins.size()), ptr);
      HypernodeID pprev = 0;
      for (HypernodeID v : pins) {
        ptr += encode_varint_to_ptr(static_cast<uint64_t>(v - pprev), ptr);
        pprev = v;
      }
    });
    // Compute totals from metadata
    hypergraph._num_pins = 0;
    hypergraph._max_edge_size = 0;
    for (HyperedgeID he = 0; he < hypergraph._num_hyperedges; ++he) {
      hypergraph._num_pins += static_cast<HypernodeID>(uniq_pins[he]);
      if (uniq_pins[he] > static_cast<uint32_t>(hypergraph._max_edge_size)) {
        hypergraph._max_edge_size = static_cast<HypernodeID>(uniq_pins[he]);
      }
    }

    // Release uniq metadata vectors we no longer need
    uniq_rep.clear(); uniq_rep.shrink_to_fit();
    uniq_weight.clear(); uniq_weight.shrink_to_fit();
    uniq_len.clear(); uniq_len.shrink_to_fit();
    uniq_pins.clear(); uniq_pins.shrink_to_fit();

  // From here on, we decode pins from the finalized incidence array using hyperedge offsets.

    // Two-pass, pre-sized encoding of incident nets to minimize peak memory
    std::vector<HyperedgeID> last_he(num_coarse_nodes, 0);
    std::vector<size_t> degree(num_coarse_nodes, 0);
    std::vector<size_t> bytes(num_coarse_nodes, 0);

    // Pass 1: compute degree and total bytes per node (including header)
    for (HyperedgeID he = 0; he < hypergraph._num_hyperedges; ++he) {
      const size_t off = hypergraph._hyperedge_offsets[he];
      const size_t end = (he + 1 < hypergraph._num_hyperedges)
                           ? hypergraph._hyperedge_offsets[he + 1]
                           : hypergraph._compressed_incidence_array.size();
      size_t pos = off;
      const uint64_t pin_cnt = decode_varint_bounded(hypergraph._compressed_incidence_array, pos, end);
      HypernodeID prev = 0;
      for (uint64_t t = 0; t < pin_cnt && pos < end; ++t) {
        const uint64_t gap = decode_varint_bounded(hypergraph._compressed_incidence_array, pos, end);
        const HypernodeID v = static_cast<HypernodeID>(prev + gap);
        prev = v;
        ++degree[v];
        const uint64_t legap = static_cast<uint64_t>(he - last_he[v]);
        bytes[v] += varint_size_u64(legap);
        last_he[v] = he;
      }
    }
    size_t total_bytes = 0;
    for (HypernodeID u = 0; u < num_coarse_nodes; ++u) {
      bytes[u] += varint_size_u64(static_cast<uint64_t>(degree[u])); // header for degree
      hypergraph._hypernode_offsets[u] = total_bytes;
      total_bytes += bytes[u];
    }

    hypergraph._compressed_incident_nets.clear();
    hypergraph._compressed_incident_nets.resize(total_bytes);

    // Initialize write pointers and write headers
    std::vector<uint8_t*> write_ptr(num_coarse_nodes, nullptr);
    uint8_t* base = hypergraph._compressed_incident_nets.data();
    for (HypernodeID u = 0; u < num_coarse_nodes; ++u) {
      uint8_t* ptr = base + hypergraph._hypernode_offsets[u];
      ptr += encode_varint_to_ptr(static_cast<uint64_t>(degree[u]), ptr);
      write_ptr[u] = ptr;
    }
    std::fill(last_he.begin(), last_he.end(), 0);

    // Pass 2: write gaps in he order
    hypergraph._total_degree = 0;
    for (HyperedgeID he = 0; he < hypergraph._num_hyperedges; ++he) {
      const size_t off = hypergraph._hyperedge_offsets[he];
      const size_t end = (he + 1 < hypergraph._num_hyperedges)
                           ? hypergraph._hyperedge_offsets[he + 1]
                           : hypergraph._compressed_incidence_array.size();
      size_t pos = off;
      const uint64_t pin_cnt = decode_varint_bounded(hypergraph._compressed_incidence_array, pos, end);
      HypernodeID prev = 0;
      for (uint64_t t = 0; t < pin_cnt && pos < end; ++t) {
        const uint64_t gapv = decode_varint_bounded(hypergraph._compressed_incidence_array, pos, end);
        const HypernodeID v = static_cast<HypernodeID>(prev + gapv);
        prev = v;
        const uint64_t legap = static_cast<uint64_t>(he - last_he[v]);
        last_he[v] = he;
        write_ptr[v] += encode_varint_to_ptr(legap, write_ptr[v]);
      }
    }
    for (HypernodeID u = 0; u < num_coarse_nodes; ++u) {
      hypergraph._total_degree += static_cast<HypernodeID>(degree[u]);
    }

  // Total node weight is sum of node weights
    HypernodeWeight tw = 0;
    for (HypernodeID u = 0; u < num_coarse_nodes; ++u) tw += hypergraph._hypernode_weights[u];
    hypergraph._total_weight = tw;

  // Release temporary buffers (none left from previous pipeline)

  // Tighten incidence buffer to actual size to reduce residual capacity
  hypergraph._compressed_incidence_array.shrink_to_fit();

    // #################### STAGE 4 ####################
    // Communities: inherit from fine graph (like static version)
    hypergraph._community_ids.assign(num_coarse_nodes, 0);
    for (HypernodeID hn = 0; hn < _num_hypernodes; ++hn) {
      const HypernodeID ch = communities[hn];
      if (ch != kInvalidHypernode) {
        hypergraph._community_ids[ch] = communityID(hn);
      }
    }

    // Fixed vertices
    if (hasFixedVertices()) {
      FixedVertexSupport<CompressedHypergraph> coarse_fixed(hypergraph.initialNumNodes(), _fixed_vertices.numBlocks());
      coarse_fixed.setHypergraph(&hypergraph);
      for (HypernodeID hn = 0; hn < _num_hypernodes; ++hn) {
        if (isFixed(hn)) {
          const HypernodeID ch = communities[hn];
          if (ch != kInvalidHypernode) {
            coarse_fixed.fixToBlock(ch, fixedVertexBlock(hn));
          }
        }
      }
      hypergraph.addFixedVertexSupport(std::move(coarse_fixed));
    }

    return hypergraph;
  }

  // ! Copy compressed hypergraph in parallel
  CompressedHypergraph CompressedHypergraph::copy(parallel_tag_t) const {
    CompressedHypergraph hypergraph;

    hypergraph._num_hypernodes = _num_hypernodes;
    hypergraph._num_removed_hypernodes = _num_removed_hypernodes;
    hypergraph._num_hyperedges = _num_hyperedges;
    hypergraph._num_removed_hyperedges = _num_removed_hyperedges;
    hypergraph._max_edge_size = _max_edge_size;
    hypergraph._num_pins = _num_pins;
    hypergraph._total_degree = _total_degree;
    hypergraph._total_weight = _total_weight;

    tbb::parallel_invoke([&] {
      hypergraph._hypernode_offsets = _hypernode_offsets;
    }, [&] {
      hypergraph._hypernode_enabled = _hypernode_enabled;
    }, [&] {
      hypergraph._compressed_incident_nets.resize(_compressed_incident_nets.size());
      memcpy(hypergraph._compressed_incident_nets.data(), _compressed_incident_nets.data(),
             sizeof(uint8_t) * _compressed_incident_nets.size());
    }, [&] {
      hypergraph._hyperedge_offsets = _hyperedge_offsets;
    }, [&] {
      hypergraph._hyperedge_enabled = _hyperedge_enabled;
    }, [&] {
      hypergraph._compressed_incidence_array.resize(_compressed_incidence_array.size());
      memcpy(hypergraph._compressed_incidence_array.data(), _compressed_incidence_array.data(),
             sizeof(uint8_t) * _compressed_incidence_array.size());
    }, [&] {
      hypergraph._community_ids = _community_ids;
    }, [&] {
      hypergraph._hypernode_weights = _hypernode_weights; // vector and map copy
    }, [&] {
      hypergraph._hyperedge_weights = _hyperedge_weights;
    }, [&] {
      hypergraph.addFixedVertexSupport(_fixed_vertices.copy());
    });
    return hypergraph;
  }

  // ! Copy compressed hypergraph sequential
  CompressedHypergraph CompressedHypergraph::copy() const {
    CompressedHypergraph hypergraph;

    hypergraph._num_hypernodes = _num_hypernodes;
    hypergraph._num_removed_hypernodes = _num_removed_hypernodes;
    hypergraph._num_hyperedges = _num_hyperedges;
    hypergraph._num_removed_hyperedges = _num_removed_hyperedges;
    hypergraph._max_edge_size = _max_edge_size;
    hypergraph._num_pins = _num_pins;
    hypergraph._total_degree = _total_degree;
    hypergraph._total_weight = _total_weight;

    hypergraph._hypernode_offsets = _hypernode_offsets;
    hypergraph._hypernode_enabled = _hypernode_enabled;

    hypergraph._compressed_incident_nets.resize(_compressed_incident_nets.size());
    // FIX: copy bytes
    memcpy(hypergraph._compressed_incident_nets.data(), _compressed_incident_nets.data(),
           sizeof(uint8_t) * _compressed_incident_nets.size());

    hypergraph._hyperedge_offsets = _hyperedge_offsets;
    hypergraph._hyperedge_enabled = _hyperedge_enabled;

    hypergraph._compressed_incidence_array.resize(_compressed_incidence_array.size());
    // FIX: copy bytes
    memcpy(hypergraph._compressed_incidence_array.data(), _compressed_incidence_array.data(),
           sizeof(uint8_t) * _compressed_incidence_array.size());

  hypergraph._community_ids = _community_ids;
  hypergraph._hypernode_weights = _hypernode_weights;
  hypergraph._hyperedge_weights = _hyperedge_weights;
    hypergraph.addFixedVertexSupport(_fixed_vertices.copy());

    return hypergraph;
  }

  void CompressedHypergraph::memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);
    parent->addChild("Hypernode Offsets", sizeof(size_t) * _hypernode_offsets.capacity());
    parent->addChild("Hypernode Enabled (bits)", (_hypernode_enabled.capacity() + 7) / 8);
    parent->addChild("Incident Nets", sizeof(uint8_t) * _compressed_incident_nets.size());
    parent->addChild("Hyperedge Offsets", sizeof(size_t) * _hyperedge_offsets.capacity());
    parent->addChild("Hyperedge Enabled (bits)", (_hyperedge_enabled.capacity() + 7) / 8);
    parent->addChild("Incidence Array", sizeof(uint8_t) * _compressed_incidence_array.size());
    parent->addChild("Communities", sizeof(PartitionID) * _community_ids.capacity());
    parent->addChild("Hypernode Weights (two-level)", _hypernode_weights.approx_bytes());
    parent->addChild("Hyperedge Weights (two-level)", _hyperedge_weights.approx_bytes());
    if ( hasFixedVertices() ) {
      parent->addChild("Fixed Vertex Support", _fixed_vertices.size_in_bytes());
    }
  }

  size_t CompressedHypergraph::memoryConsumptionKB() const {
    size_t total = 0;
    total += sizeof(size_t) * _hypernode_offsets.capacity();
    // BitVector packs bits; approximate bytes as ceil(capacity_in_bits/8)
    total += ((_hypernode_enabled.capacity() + 7) / 8);
    total += sizeof(uint8_t) * _compressed_incident_nets.size();
    total += sizeof(size_t) * _hyperedge_offsets.capacity();
    total += ((_hyperedge_enabled.capacity() + 7) / 8);
    total += sizeof(uint8_t) * _compressed_incidence_array.size();
    total += _hypernode_weights.approx_bytes();
    total += _hyperedge_weights.approx_bytes();
    total += sizeof(PartitionID) * _community_ids.capacity();
    if (_fixed_vertices.hasFixedVertices()) {
        total += _fixed_vertices.size_in_bytes();
    }
    return total / 1024;
  }

  // ! Computes the total node weight of the hypergraph
  void CompressedHypergraph::computeAndSetTotalNodeWeight(parallel_tag_t) {
    _total_weight = tbb::parallel_reduce(tbb::blocked_range<HypernodeID>(ID(0), _num_hypernodes), 0,
                                         [this](const tbb::blocked_range<HypernodeID>& range, HypernodeWeight init) {
                                           HypernodeWeight weight = init;
                                           for (HypernodeID hn = range.begin(); hn < range.end(); ++hn) {
                                             if (nodeIsEnabled(hn)) {
                                               weight += this->nodeWeight(hn);
                                             }
                                           }
                                           return weight;
                                         }, std::plus<>());
  }

} // namespace

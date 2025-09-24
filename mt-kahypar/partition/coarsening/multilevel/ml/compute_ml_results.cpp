/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2025 Nikolai Maas <nikolai.maas@kit.edu>
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

#include "compute_ml_results.h"

#include <memory>

#include <tbb/parallel_invoke.h>
#include <tbb/enumerable_thread_specific.h>

#include "mt-kahypar/partition/coarsening/multilevel/ml/compute_features.h"
#include "mt-kahypar/partition/coarsening/multilevel/ml/feature_definitions.h"
#include "mt-kahypar/partition/coarsening/multilevel/ml/mlp.h"

#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

#define is_aligned(POINTER, BYTE_COUNT) \
    (((uintptr_t)(const void *)(POINTER)) % (BYTE_COUNT) == 0)

constexpr size_t EDGE_CHUNK_SIZE = 1024;

struct InOutPredictionBuffer {
  float input_buffer alignas(32) [(EDGE_CHUNK_SIZE + 4) * 48];
  float output_buffer alignas(32) [(EDGE_CHUNK_SIZE + 4) * 48];
};

using BufferPtr = std::unique_ptr<InOutPredictionBuffer>;


template<typename T, size_t DEPTH>
struct MapFeaturesForEdgeImpl;

template<typename MAPPER, typename... TAIL, size_t DEPTH>
struct MapFeaturesForEdgeImpl<kahypar::meta::Typelist<MAPPER, TAIL...>, DEPTH> {
  static void map_features(const GlobalFeatures& global, const N1Features& n1_features_0, const N1Features& n1_features_1, const EdgeFeatures& edge_features, float* output) {
    *output = MAPPER::map(global, n1_features_0, n1_features_1, edge_features);
    ++output;
    MapFeaturesForEdgeImpl<kahypar::meta::Typelist<TAIL...>, DEPTH + 1>::map_features(global, n1_features_0, n1_features_1, edge_features, output);
  }
};

template<size_t DEPTH>
struct MapFeaturesForEdgeImpl<kahypar::meta::Typelist<>, DEPTH> {
  static void map_features(const GlobalFeatures&, const N1Features&, const N1Features&, const EdgeFeatures&, float*) {
    static_assert(DEPTH == 32);
  }
};

using MapFeaturesForEdge = MapFeaturesForEdgeImpl<features::OrderedFeatures, 0>;


MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
void predictEdgesImpl(const ds::StaticGraph& graph, const Context& context,
                      const GlobalFeatures& global, const ds::Array<N1Features>& n1_features, bool skip_comm_1,
                      HyperedgeID first_edge, HyperedgeID n,
                      float *params_row, InOutPredictionBuffer& buffer, vec<EdgeMetadata>& metadata) {
  float* data = buffer.input_buffer;
  for (HyperedgeID he = first_edge; he < first_edge + n; ++he) {
    ASSERT(graph.edgeIsEnabled(he));

    HypernodeID u = graph.edgeSource(he);
    HypernodeID v = graph.edgeTarget(he);
    const N1Features* u_f = &n1_features[u];
    const N1Features* v_f = &n1_features[v];

    // TODO: don't predict edges twice??
    if (u_f->degree < v_f->degree || (u_f->degree == v_f->degree && u < v)) {
      std::swap(u, v);
      std::swap(u_f, v_f);
    }
    EdgeFeatures edge = computeEdgeFeatures(graph, context, u, *u_f, v, *v_f, skip_comm_1);
    MapFeaturesForEdge::map_features(global, *u_f, *v_f, edge, data);

    static_assert(MEANS.size() == STDEVS.size() && MEANS.size() == 32);
    for (size_t i = 0; i < MEANS.size(); ++i) {
      float val = data[i];
      data[i] = (val - MEANS[i]) / STDEVS[i];
    }

    data += 32;
  }

  float* in = buffer.input_buffer;
  float* out = buffer.output_buffer;
  ASSERT(is_aligned(in, 32) && is_aligned(out, 32));
  predict(in, params_row, n, out);

  for (HyperedgeID offset = 0; offset < n; ++offset) {
    ASSERT(out[offset] >= 0 && out[offset] <= 1, V(out[offset]));
    metadata[first_edge + offset] = out[offset];
  }
}

void predictEdgesStaticSize(const ds::StaticGraph& graph, const Context& context,
                            const GlobalFeatures& global, const ds::Array<N1Features>& n1_features, bool skip_comm_1,
                            HyperedgeID first_edge,
                            float *params_row, InOutPredictionBuffer& buffer, vec<EdgeMetadata>& metadata) {
  predictEdgesImpl(graph, context, global, n1_features, skip_comm_1, first_edge, EDGE_CHUNK_SIZE, params_row, buffer, metadata);
}

void predictEdgesDynamicSize(const ds::StaticGraph& graph, const Context& context,
                             const GlobalFeatures& global, const ds::Array<N1Features>& n1_features, bool skip_comm_1,
                             HyperedgeID first_edge, HyperedgeID n,
                             float *params_row, InOutPredictionBuffer& buffer, vec<EdgeMetadata>& metadata) {
  predictEdgesImpl(graph, context, global, n1_features, skip_comm_1, first_edge, n, params_row, buffer, metadata);
}

void computeEdgeMetadataFromModel(const ds::StaticGraph& graph, const Context& context, vec<EdgeMetadata>& metadata) {
  GlobalFeatures global_features;
  ds::Array<N1Features> n1_features;
  bool skip_comm_1;

  tbb::parallel_invoke([&] {
      std::tie(global_features, n1_features, skip_comm_1) = computeFeatures(graph, context);
    }, [&] {
      metadata.resize(graph.initialNumEdges());
    }
  );

  tbb::enumerable_thread_specific<BufferPtr> in_out_buffers([]{
    return BufferPtr(new InOutPredictionBuffer());
  });
  float* params_row = precompute_params();

  utils::Timer& timer = utils::Utilities::instance().getTimer(context.utility_id);
  timer.start_timer("inference", "Model Inference");

  size_t n_chunks = (graph.initialNumEdges() + EDGE_CHUNK_SIZE - 1) / EDGE_CHUNK_SIZE;
  tbb::parallel_for(UL(0), n_chunks, [&](size_t chunk_id) {
    HyperedgeID first_edge = chunk_id * EDGE_CHUNK_SIZE;
    InOutPredictionBuffer& buffer = *in_out_buffers.local();
    if (first_edge + EDGE_CHUNK_SIZE < graph.initialNumEdges()) {
      predictEdgesStaticSize(graph, context, global_features, n1_features, skip_comm_1, first_edge, params_row, buffer, metadata);
    } else {
      size_t n = graph.initialNumEdges() - first_edge;
      predictEdgesDynamicSize(graph, context, global_features, n1_features, skip_comm_1, first_edge, n, params_row, buffer, metadata);
    }
  });

  timer.stop_timer("inference");
  free_params(params_row);
}

}  // namespace mt_kahypar

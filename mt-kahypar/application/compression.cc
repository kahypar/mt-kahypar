#include <chrono>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
namespace fs = std::filesystem;

#include <iostream>
#include <algorithm>
#include <cstdlib>
#ifdef __linux__
#include <malloc.h>
#endif

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/io/presets.h"
#include "mt-kahypar/parallel/tbb_initializer.h"
#include "mt-kahypar/partition/conversion.h"
#include "mt-kahypar/partition/mapping/target_graph.h"
#include "mt-kahypar/partition/partitioner_facade.h"
#include "mt-kahypar/partition/registries/register_memory_pool.h"
#include "mt-kahypar/partition/registries/registry.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/delete.h"
#include "mt-kahypar/utils/exception.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/utilities.h"

#include "mt-kahypar/utils/memory_tree.h"

#include "mt-kahypar/datastructures/compressed_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"

using namespace mt_kahypar;

// Simple partition quality estimator for comparisons between static and compressed outputs.
// We mimic common objectives without depending on metrics.h:
// - cut: sum edge weight for edges with connectivity > 1
// - km1: sum (connectivity - 1) * edge weight
// - soed: sum connectivity * edge weight for edges with connectivity > 1
// - steiner_tree: approximate via (connectivity - 1) * edge weight
template <typename PartitionedHypergraphT>
static inline HyperedgeWeight simple_partition_quality(const PartitionedHypergraphT& phg,
                                                       const Objective objective) {
    HyperedgeWeight total = 0;
    switch (objective) {
        case Objective::cut: {
            for (const HyperedgeID he : phg.edges()) {
                const PartitionID c = phg.connectivity(he);
                if (c > 1) total += phg.edgeWeight(he);
            }
            break;
        }
        case Objective::km1: {
            for (const HyperedgeID he : phg.edges()) {
                const PartitionID c = phg.connectivity(he);
                if (c > 1) total += static_cast<HyperedgeWeight>((c - 1)) * phg.edgeWeight(he);
            }
            break;
        }
        case Objective::soed: {
            for (const HyperedgeID he : phg.edges()) {
                const PartitionID c = phg.connectivity(he);
                if (c > 1) total += static_cast<HyperedgeWeight>(c) * phg.edgeWeight(he);
            }
            break;
        }
        case Objective::steiner_tree: {
            // Cheap approximation without target graph: treat like km1
            for (const HyperedgeID he : phg.edges()) {
                const PartitionID c = phg.connectivity(he);
                if (c > 1) total += static_cast<HyperedgeWeight>((c - 1)) * phg.edgeWeight(he);
            }
            break;
        }
        default: break;
    }
    // Graphs store each edge twice; hypergraphs do not. Mirror metrics behavior.
    if constexpr (PartitionedHypergraphT::is_graph) {
        total /= 2;
    }
    return total;
}

static inline HyperedgeWeight compute_partition_quality_estimate(
    const mt_kahypar_partitioned_hypergraph_t& phg_handle, const Context& context) {
    switch (phg_handle.type) {
        case MULTILEVEL_HYPERGRAPH_PARTITIONING: {
            const StaticPartitionedHypergraph& phg =
                utils::cast<StaticPartitionedHypergraph>(phg_handle);
            return simple_partition_quality(phg, context.partition.objective);
        }
        case COMPRESSED_MULTILEVEL_HYPERGRAPH_PARTITIONING: {
            const CompressedPartitionedHypergraph& phg =
                utils::cast<CompressedPartitionedHypergraph>(phg_handle);
            return simple_partition_quality(phg, context.partition.objective);
        }
        #ifdef KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES
        case N_LEVEL_HYPERGRAPH_PARTITIONING: {
            const DynamicPartitionedHypergraph& phg =
                utils::cast<DynamicPartitionedHypergraph>(phg_handle);
            return simple_partition_quality(phg, context.partition.objective);
        }
        #endif
        #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
        case MULTILEVEL_GRAPH_PARTITIONING: {
            const StaticPartitionedGraph& pg = utils::cast<StaticPartitionedGraph>(phg_handle);
            return simple_partition_quality(pg, context.partition.objective);
        }
        case N_LEVEL_GRAPH_PARTITIONING: {
            const DynamicPartitionedGraph& pg = utils::cast<DynamicPartitionedGraph>(phg_handle);
            return simple_partition_quality(pg, context.partition.objective);
        }
        #endif
        #ifdef KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES
        case LARGE_K_PARTITIONING:
        #endif
        case NULLPTR_PARTITION:
        default:
        break;
    }
    return 0;
}

// ---- Memory usage helpers (Linux) ----
static inline size_t read_proc_status_kb(const char* key) {
    std::ifstream in("/proc/self/status");
    if (!in) return 0;
    std::string line;
    const std::string needle = std::string(key) + ":";
    while (std::getline(in, line)) {
        if (line.compare(0, needle.size(), needle) == 0) {
            // Format: Key:\t<value> kB
            size_t kb = 0;
            for (size_t i = needle.size(); i < line.size(); ++i) {
                if (line[i] >= '0' && line[i] <= '9') {
                    kb = std::strtoull(line.c_str() + i, nullptr, 10);
                    return kb;
                }
            }
        }
    }
    return 0;
}

static inline size_t getCurrentRSSKB() { return read_proc_status_kb("VmRSS"); }
static inline size_t getPeakRSSKB()    { return read_proc_status_kb("VmHWM"); }

// No sampler needed: we use VmHWM (MaxRSS) checkpoints per phase and process isolation between files/variants.

// Write a fresh header (truncate file). Call this once per full benchmark run (parent process).
void write_csv_header(const std::string& csv_path) {
    std::ofstream outfile(csv_path, std::ios::trunc);
    outfile << "FileName,NodeCount,HyperedgeCount,IOTimeMS,PartitionTimeMS,"
               "IOPeakRSSKB,PartitionPeakRSSKB,IOAbsPeakRSSKB,PartitionAbsPeakRSSKB,"
               "MemoryUsage,PartitionQuality,Compressed"
            << std::endl;
}

// Ensure header exists without truncation (used by per-file child runs).
static inline void ensure_csv_header(const std::string& csv_path) {
    std::ifstream in(csv_path);
    bool need_header = true;
    if (in) {
        std::string first;
        if (std::getline(in, first)) {
            if (first.rfind("FileName,NodeCount,HyperedgeCount,IOTimeMS,PartitionTimeMS,", 0) == 0) {
                need_header = false;
            }
        }
    }
    if (need_header) {
        std::ofstream outfile(csv_path, std::ios::app);
        outfile << "FileName,NodeCount,HyperedgeCount,IOTimeMS,PartitionTimeMS,"
                   "IOPeakRSSKB,PartitionPeakRSSKB,IOAbsPeakRSSKB,PartitionAbsPeakRSSKB,"
                   "MemoryUsage,PartitionQuality,Compressed"
                << std::endl;
    }
}

void append_csv_line(const std::string& csv_path, const std::string& filename,
                     size_t num_nodes, size_t num_hyperedges, size_t io_time,
                     size_t partition_time,
                     size_t io_peak_kb, size_t part_peak_kb,
                     size_t io_abs_peak_kb, size_t part_abs_peak_kb,
                     size_t memory_usage,
                     HyperedgeWeight partition_quality,
                     bool compressed) {
    std::ofstream outfile(csv_path, std::ios::app);
    outfile << filename << "," << num_nodes << "," << num_hyperedges << ","
            << io_time << "," << partition_time << ","
            << io_peak_kb << "," << part_peak_kb << ","
            << io_abs_peak_kb << "," << part_abs_peak_kb << ","
            << memory_usage << ","
            << partition_quality << ","
            << (compressed ? "1" : "0") << std::endl;
}

void partition_graph(const std::string& filename, const Context& base_context,
                     const std::string& OUTPUT_PATH, bool compressed = false) {
    // Make a fresh, per-run copy of the context
    Context context = base_context;

    // Ensure Context invariant before sanityCheck
    if (!context.partition.use_individual_part_weights &&
        !context.partition.max_part_weights.empty()) {
      context.partition.max_part_weights.clear();
    }

    context.partition.graph_filename = filename;
    context.partition.instance_type = compressed
                                          ? InstanceType::compressed_hypergraph
                                          : InstanceType::hypergraph;
    context.partition.partition_type = to_partition_c_type(
        context.partition.preset_type, context.partition.instance_type);

    context.utility_id =
        utils::Utilities::instance().registerNewUtilityObjects();
    if (context.partition.verbose_output) io::printBanner();

    utils::Randomize::instance().setSeed(context.partition.seed);
    if (context.shared_memory.use_localized_random_shuffle) {
        utils::Randomize::instance().enableLocalizedParallelShuffle(
            context.shared_memory.shuffle_block_size);
    }

#ifndef KAHYPAR_DISABLE_HWLOC
    size_t num_available_cpus = HardwareTopology::instance().num_cpus();
    if (num_available_cpus < context.shared_memory.num_threads) {
        WARNING("There are currently only"
                << num_available_cpus << "cpus available."
                << "Setting number of threads from"
                << context.shared_memory.num_threads << "to"
                << num_available_cpus);
        context.shared_memory.num_threads = num_available_cpus;
    }
#endif

    // Initialize TBB task arenas on numa nodes
    TBBInitializer::instance(context.shared_memory.num_threads);

#ifndef KAHYPAR_DISABLE_HWLOC
    // We set the membind policy to interleaved allocations in order to
    // distribute allocations evenly across NUMA nodes
    hwloc_cpuset_t cpuset = TBBInitializer::instance().used_cpuset();
    parallel::HardwareTopology<>::instance()
        .activate_interleaved_membind_policy(cpuset);
    hwloc_bitmap_free(cpuset);
#endif

		HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    // Read Hypergraph
    utils::Timer& timer =
        utils::Utilities::instance().getTimer(context.utility_id);
    // Memory baseline before I/O
    const size_t io_start_rss_kb = getCurrentRSSKB();
    timer.start_timer("io_hypergraph", "I/O Hypergraph");

    mt_kahypar_hypergraph_t hypergraph = io::readInputFile(
        context.partition.graph_filename, context.partition.preset_type,
        context.partition.instance_type, context.partition.file_format,
        context.preprocessing.stable_construction_of_incident_edges);

    timer.stop_timer("io_hypergraph");
    // Memory after I/O: use MaxRSS (VmHWM) checkpoint
    const size_t io_abs_peak_kb = getPeakRSSKB();
    const size_t io_end_rss_kb = getCurrentRSSKB();
    const size_t io_rss_growth = (io_end_rss_kb >= io_start_rss_kb) ? (io_end_rss_kb - io_start_rss_kb) : 0;
    size_t io_peak_kb = (io_abs_peak_kb > io_start_rss_kb) ? (io_abs_peak_kb - io_start_rss_kb) : io_rss_growth;
		HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
		size_t io_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    // Read Target Graph
    std::unique_ptr<TargetGraph> target_graph;
    if (context.partition.objective == Objective::steiner_tree) {
        if (context.mapping.target_graph_file != "") {
            target_graph = std::make_unique<TargetGraph>(
                io::readInputFile<ds::StaticGraph>(
                    context.mapping.target_graph_file, FileFormat::Metis,
                    true));
        } else {
            throw InvalidInputException(
                "No target graph file specified (use -g <file> or "
                "--target-graph-file=<file>)!");
        }
    }

    if (context.partition.fixed_vertex_filename != "") {
        timer.start_timer("read_fixed_vertices", "Read Fixed Vertex File");
        io::addFixedVerticesFromFile(hypergraph,
                                     context.partition.fixed_vertex_filename,
                                     context.partition.k);
        timer.stop_timer("read_fixed_vertices");
    }

    // Initialize Memory Pool and Algorithm/Policy Registries
    register_memory_pool(hypergraph, context);
    register_algorithms_and_policies();

    // Partition Hypergraph
    start = std::chrono::high_resolution_clock::now();
    const size_t part_start_rss_kb = getCurrentRSSKB();
    mt_kahypar_partitioned_hypergraph_t partitioned_hypergraph =
        PartitionerFacade::partition(hypergraph, context, target_graph.get());
    end = std::chrono::high_resolution_clock::now();
    const size_t part_abs_peak_kb = getPeakRSSKB();
    const size_t part_end_rss_kb = getCurrentRSSKB();

		size_t partition_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    // Compute a simple partition quality estimate (no dependency on metrics.h)
    const HyperedgeWeight partition_quality_estimate =
        compute_partition_quality_estimate(partitioned_hypergraph, context);

    // Print Stats
    std::chrono::duration<double> elapsed_seconds(end - start);
    PartitionerFacade::printPartitioningResults(partitioned_hypergraph, context,
                                                elapsed_seconds);
    if (context.partition.sp_process_output) {
        std::cout << PartitionerFacade::serializeResultLine(
                         partitioned_hypergraph, context, elapsed_seconds)
                  << std::endl;
    }
    if (context.partition.csv_output) {
        std::cout << PartitionerFacade::serializeCSV(partitioned_hypergraph,
                                                     context, elapsed_seconds)
                  << std::endl;
    }
    if (context.partition.write_partition_file) {
        PartitionerFacade::writePartitionFile(
            partitioned_hypergraph, context.partition.graph_partition_filename);
    }
    // Optional one-shot memory tree dump per variant, controlled by env var MTK_MEMTREE_ONCE
    auto env_memtree = std::getenv("MTK_MEMTREE_ONCE");

    if (compressed) {
        const ds::CompressedHypergraph& hg =
            utils::cast<ds::CompressedHypergraph>(hypergraph);
        if (env_memtree && *env_memtree != '\0') {
            static bool dumped_compressed_once = false;
            if (!dumped_compressed_once) {
                utils::MemoryTreeNode root("CompressedHypergraph");
                hg.memoryConsumption(&root);
                root.finalize();
                std::ofstream mt_file("./__out/memory_tree_compressed.txt", std::ios::out | std::ios::trunc);
                if (mt_file) { mt_file << root; }
                dumped_compressed_once = true;
            }
        }
                append_csv_line(OUTPUT_PATH, filename, hg.initialNumNodes(),
                                                hg.initialNumEdges(), io_time, partition_time,
                                                                        io_peak_kb,
                                                                        (part_abs_peak_kb > part_start_rss_kb ? (part_abs_peak_kb - part_start_rss_kb) :
                                                                            (part_end_rss_kb >= part_start_rss_kb ? part_end_rss_kb - part_start_rss_kb : 0UL)
                                                                        ),
                                                                        io_abs_peak_kb,
                                                                        part_abs_peak_kb,
                                                hg.memoryConsumptionKB(),
                        partition_quality_estimate,
                        compressed);
    } else {
        const ds::StaticHypergraph& hg =
            utils::cast<ds::StaticHypergraph>(hypergraph);
        if (env_memtree && *env_memtree != '\0') {
            static bool dumped_static_once = false;
            if (!dumped_static_once) {
                utils::MemoryTreeNode root("StaticHypergraph");
                hg.memoryConsumption(&root);
                root.finalize();
                std::ofstream mt_file("./__out/memory_tree_static.txt", std::ios::out | std::ios::trunc);
                if (mt_file) { mt_file << root; }
                dumped_static_once = true;
            }
        }
                append_csv_line(OUTPUT_PATH, filename, hg.initialNumNodes(),
                                                hg.initialNumEdges(), io_time, partition_time,
                                                                        io_peak_kb,
                                                                        (part_abs_peak_kb > part_start_rss_kb ? (part_abs_peak_kb - part_start_rss_kb) :
                                                                            (part_end_rss_kb >= part_start_rss_kb ? part_end_rss_kb - part_start_rss_kb : 0UL)
                                                                        ),
                                                                        io_abs_peak_kb,
                                                                        part_abs_peak_kb,
                                                hg.memoryConsumptionKB(),
                        partition_quality_estimate,
                        compressed);
    }

        parallel::MemoryPool::instance().free_memory_chunks();
        // Optionally ask glibc to return freed memory to the OS to avoid baseline inflation across files
        if (const char* trim = std::getenv("MTK_MALLOC_TRIM"); trim && *trim == '1') {
            #ifdef __linux__
            malloc_trim(0);
            #endif
        }
    TBBInitializer::instance().terminate();

    utils::delete_hypergraph(hypergraph);
    utils::delete_partitioned_hypergraph(partitioned_hypergraph);
}

void print_progress(size_t current, size_t total) {
    const int bar_width = 50;  // Width of the progress bar
    float progress = static_cast<float>(current) / total;

    std::cout << "\r[";
    int pos = static_cast<int>(bar_width * progress);
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << current << "/" << total << " files" << std::flush;
}

int main(int argc, char* argv[]) {
    Context base(false);
    processCommandLineInput(base, argc, argv, nullptr);
    base.partition.verbose_output = false;

    fs::create_directories("./__out");

    const std::string OUTPUT_PATH = "./__out/benchmark.csv";

    // Child mode: run a single file (used for isolation via environment variables)
    if (const char* single = std::getenv("MTK_SINGLE_FILE")) {
        const char* c = std::getenv("MTK_COMPRESSED");
        const bool compressed = (c && *c == '1');
        ensure_csv_header(OUTPUT_PATH);
        partition_graph(single, base, OUTPUT_PATH, compressed);
        return 0;
    }

    // Parent process: reset CSV for a fresh run
    write_csv_header(OUTPUT_PATH);

    const std::string directory = "./_graphs/benchmark_set_d";

    // Collect all regular files first
    std::vector<fs::path> files;
    for (const auto& entry : fs::directory_iterator(directory)) {
        if (!fs::is_regular_file(entry)) continue;
        files.push_back(entry.path());
    }

    size_t total = files.size();
    if (total == 0) {
        std::cout << "No files found in " << directory << std::endl;
        return 1;
    }

    // Always run isolated: spawn a fresh process per file and variant to avoid carry-over
    {
        // Reconstruct base command (program + original args)
        std::string base_cmd;
        {
            // Quote arguments to be safe for spaces
            auto q = [](const std::string& s){
                std::string r;
                r.reserve(s.size() + 2);
                r.push_back('\'');
                for (char ch : s) {
                    if (ch == '\'') { r += "'\\''"; }
                    else { r.push_back(ch); }
                }
                r.push_back('\'');
                return r;
            };
            base_cmd = q(argv[0]);
            for (int i = 1; i < argc; ++i) {
                base_cmd += " ";
                base_cmd += q(argv[i]);
            }
        }

        size_t count = 0;
        std::cout << std::endl << "Partitioning Uncompressed Graphs (isolated)" << std::endl;
        print_progress(count, total);
        for (const auto& path : files) {
            std::string cmd = std::string("MTK_SINGLE_FILE=") + path.string() +
                              " MTK_COMPRESSED=0 MTK_MALLOC_TRIM=1 " + base_cmd + " > /dev/null";
            int rc = std::system(cmd.c_str());
            (void)rc;
            print_progress(++count, total);
        }
        std::cout << std::endl << std::endl <<  "Partitioning Compressed Graphs (isolated)" << std::endl;
        count = 0;
        print_progress(count, total);
        for (const auto& path : files) {
            std::string cmd = std::string("MTK_SINGLE_FILE=") + path.string() +
                              " MTK_COMPRESSED=1 MTK_MALLOC_TRIM=1 " + base_cmd + " > /dev/null";
            int rc = std::system(cmd.c_str());
            (void)rc;
            print_progress(++count, total);
        }
        std::cout << std::endl << "Processing complete." << std::endl;
        return 0;
    }
}
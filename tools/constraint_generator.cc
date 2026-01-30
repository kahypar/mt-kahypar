#include <boost/program_options.hpp>

#include <optional>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <random>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"

using namespace mt_kahypar;
namespace po = boost::program_options;
namespace fs = std::filesystem;

const std::string CONSTRAINT_FILE_EXTENSION = ".constraints.txt";

std::optional<fs::path> findFileWithPrefix(const fs::path& path, const std::string& prefix) {
    if (fs::is_regular_file(path)) {
        auto name = path.filename().string();
        if (name.rfind(prefix, 0) == 0) {  // starts with prefix
            return path;
        }
    }
    if (fs::is_directory(path)) {
        for (const auto& entry : fs::directory_iterator(path)) {
            if (entry.is_regular_file()) {
                auto filename = entry.path().filename().string();
                if (filename.rfind(prefix, 0) == 0) {
                    return entry.path();
                }
            }
        }
    }
    return std::nullopt;
}

void print_constraints(fs::path constraints_path, 
                        vec<std::pair<HypernodeID, HypernodeID>>& constraint_list, 
                        HypernodeID num_constraints) {
    std::mt19937 gen(12345);
    std::shuffle(constraint_list.begin(), constraint_list.end(), gen);

    std::ofstream out_stream(constraints_path);
    for (HypernodeID i = 0; i < num_constraints; i++) {
        std::pair<HypernodeID, HypernodeID> constraint = constraint_list[i];
        out_stream << constraint.first << " " << constraint.second << std::endl;
    }
    out_stream.close();
}

HypernodeID pick_random(const HypernodeID limit) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, limit);
    return dist(gen);
}

HypernodeID generate_constraints_from_hg(const fs::path hg_path,
                                            const fs::path constraints_path,
                                            const float constraints_percentage,
                                            const HypernodeID max_constraints_per_node) {
    // Read Hypergraph
    HyperedgeID num_edges;
    HypernodeID num_nodes;
    HyperedgeID num_removed_single_pin_hyperedges;
    io::HyperedgeVector hyperedges;
    vec<HyperedgeWeight> hyperedges_weight;
    vec<HypernodeWeight> hypernodes_weight;

    io::readHypergraphFile(hg_path.string(), num_edges, num_nodes, num_removed_single_pin_hyperedges,
                         hyperedges, hyperedges_weight, hypernodes_weight);
    ALWAYS_ASSERT(hyperedges.size() == num_edges);
    ALWAYS_ASSERT(num_removed_single_pin_hyperedges == 0);

    HypernodeID num_constraints = num_nodes * constraints_percentage;
    HypernodeID constraint_count = 0;
    std::unordered_map<HypernodeID, HypernodeID> constraints_per_node;

    std::ofstream out_stream(constraints_path);
    while(constraint_count < num_constraints) {
        HypernodeID node = pick_random(num_nodes - 1);
        HypernodeID node_constraints = pick_random(max_constraints_per_node);
        HypernodeID count_node_constraints = 0;
        while (count_node_constraints < node_constraints && constraint_count < num_constraints) {
            HypernodeID other_node = pick_random(num_nodes - 1);
            if (constraints_per_node[node] >= max_constraints_per_node) break;
            if (other_node == node || constraints_per_node[other_node] >= max_constraints_per_node) continue;

            count_node_constraints++;
            constraint_count++;
            constraints_per_node[node]++;
            constraints_per_node[other_node]++;
            out_stream << node << " " << other_node << std::endl;
        }
    }
    out_stream.close();
    return constraint_count;
}

HypernodeID generate_constraints_from_partitioned_hg(const fs::path hg_path, 
                                                        const fs::path part_hg_path, 
                                                        const fs::path constraints_path, 
                                                        const float constraints_percentage, 
                                                        const HypernodeID max_constraints_per_node) {
    // Read Hypergraph header and partitioned Hypergraph
    HyperedgeID num_edges;
    HypernodeID num_nodes;
    std::vector<PartitionID> partitions;
    io::onlyReadHGRHeader(hg_path.string(), num_edges, num_nodes);
    io::readPartitionFile(part_hg_path.string(), num_nodes, partitions);

    HypernodeID num_constraints = num_nodes * constraints_percentage;

    vec<std::pair<HypernodeID, HypernodeID>> constraint_list;
    constraint_list.reserve(num_nodes * max_constraints_per_node);
    HypernodeID constraint_count = 0;
    std::unordered_map<HypernodeID, HypernodeID> constraints_per_node;
    for (HypernodeID node = 0; node < num_nodes; node++) {
        PartitionID node_partition = partitions[node];
        for (HypernodeID other_node = node + 1; other_node < num_nodes; other_node++) {
            if (constraints_per_node[node] >= max_constraints_per_node) break;
            if (constraints_per_node[other_node] >= max_constraints_per_node) continue;

            if (node_partition != partitions[other_node]) {
                constraint_count++;
                constraints_per_node[node]++;
                constraints_per_node[other_node]++;
                constraint_list.push_back(std::make_pair(node, other_node));
            }
        }
        constraints_per_node.erase(node);
    }
    HypernodeID actual_num_constraints = std::min(num_constraints, constraint_count);
    print_constraints(constraints_path, constraint_list, actual_num_constraints);
    return actual_num_constraints;
}

int main(int argc, char* argv[]) {
    fs::path hypergraph_path;
    vec<fs::path> hypergraph_files;
    std::optional<fs::path> partitioned_hypergraph_path;
    fs::path constraint_dir;
    HypernodeID k;
    float num_constraints_percentage;
    HypernodeID max_constraints_per_node;

    po::options_description options("Options");
    options.add_options()
        ("hypergraph,h",
        po::value<fs::path>(&hypergraph_path)->value_name("<path>")->required(),
        "Hypergraph Filename or directory")
        ("constraint,c",
        po::value<fs::path>(&constraint_dir)->value_name("<path>")->required(),
        "Constraint directory")
        ("blocks,k",
        po::value<HypernodeID>(&k)->value_name("<int>")->required(),
        "Number of blocks")
        ("num-constraints,n",
        po::value<float>(&num_constraints_percentage)->value_name("<float>")->default_value(1.0),
        "Number of constraints (optional)")
        ("part-hypergraph,p",
        po::value<fs::path>()->notifier([&](const fs::path& path) {
            partitioned_hypergraph_path = path;
        })->value_name("<path>"),
        "Partitioned-Hypergraph Filename or directory (optional)")
    ;
    po::variables_map cmd_vm;
    po::store(po::parse_command_line(argc,argv, options), cmd_vm);
    po::notify(cmd_vm);

    max_constraints_per_node = k - 1;

    if (!fs::is_directory(constraint_dir)) {
        throw fs::filesystem_error(
            "constraint path is not a directory",
            constraint_dir,
            std::make_error_code(std::errc::not_a_directory)
        );
    }

    if (fs::is_regular_file(hypergraph_path) && hypergraph_path.extension() == ".hgr") {
        hypergraph_files.push_back(hypergraph_path);
    } else if (fs::is_directory(hypergraph_path)) {
        for (const auto& file : fs::directory_iterator(hypergraph_path)) {
            if (fs::is_regular_file(file.path()) && file.path().extension() == ".hgr") {
                hypergraph_files.push_back(file.path());
            }
        }
    } else {
        throw fs::filesystem_error(
            "hypergraph file doesnt exist",
            hypergraph_path,
            std::make_error_code(std::errc::not_a_directory)
        );
    }

    for (const auto& hg_file : hypergraph_files) {
        fs::path constraint_file = constraint_dir / (hg_file.filename().string() + "." + std::to_string(k) + CONSTRAINT_FILE_EXTENSION);
        HypernodeID generated_constraints;
        if (partitioned_hypergraph_path) {
            auto part_hg_file = findFileWithPrefix(partitioned_hypergraph_path.value(), hg_file.filename().string() + ".part" + std::to_string(k));
            if (!part_hg_file){
                throw fs::filesystem_error(
                    "no matching partitioned hypergraph file found in ",
                    partitioned_hypergraph_path.value(),
                    std::make_error_code(std::errc::not_a_directory)
                );
            }
            generated_constraints = generate_constraints_from_partitioned_hg(hg_file, part_hg_file.value(), constraint_file, num_constraints_percentage, max_constraints_per_node);
        } else {
            generated_constraints = generate_constraints_from_hg(hg_file, constraint_file, num_constraints_percentage, max_constraints_per_node);
        }
        LOG << "";
        LOG << "Generated " << generated_constraints << "constraints for hg:" << hg_file.filename();
    }
    LOG << "";
    return 0;
}
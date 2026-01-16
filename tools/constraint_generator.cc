#include <boost/program_options.hpp>

#include <optional>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <random>
#include <unordered_set>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"

using namespace mt_kahypar;
namespace po = boost::program_options;
namespace fs = std::filesystem;

constexpr double DEFAULT_CONSTRAINT_FRACTION = 1.0;
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

void logNumberOfDifferentNodes(fs::path constraints_path, 
                        vec<std::pair<HypernodeID, HypernodeID>>& constraint_list, 
                        HypernodeID num_constraints) {
    std::mt19937 gen(12345);
    std::shuffle(constraint_list.begin(), constraint_list.end(), gen);
    std::unordered_set<HypernodeID> set;
    set.reserve(constraint_list.size());

    for (HypernodeID i = 0; i < num_constraints; i++) {
        std::pair<HypernodeID, HypernodeID> constraint = constraint_list[i];
        set.insert(constraint.first);
        set.insert(constraint.second);
    }
    LOG << "Number of different Nodes:"<<set.size();
}

HypernodeID get_constraints(vec<vec<NodeID>>& adjacency, 
                            const io::HyperedgeVector& hyperedges, 
                            const HypernodeID max_constraints_per_node, 
                            const HypernodeID num_constraints) {
    HypernodeID constraint_count = 0;
    std::unordered_map<HypernodeID, HypernodeID> constraints_per_node;
    for (io::Hyperedge edge : hyperedges) {
        for (HyperedgeID i = 0; i < edge.size(); i++) {
            HypernodeID node = edge[i];
            for (HyperedgeID j = i + 1; j < edge.size(); j++) {
                HypernodeID other_node = edge[j];
                if (constraint_count >= num_constraints){
                    return constraint_count;
                }
                if (constraints_per_node[node] >= max_constraints_per_node) break;
                if (constraints_per_node[other_node] >= max_constraints_per_node) continue;
                
                auto& list = adjacency[node];
                // if constraint already exists
                if (std::find(list.begin(), list.end(), other_node) != list.end()) continue;
                auto& other_list = adjacency[other_node];
                if (std::find(other_list.begin(), other_list.end(), node) != other_list.end()) continue;
                
                constraint_count++;
                constraints_per_node[node]++;
                constraints_per_node[other_node]++;
                adjacency[node].push_back(other_node);
            }
        }
    }
    return constraint_count;
}

HypernodeID generate_constraints_from_hg(const fs::path hg_path,
                                            const fs::path constraints_path,
                                            HypernodeID num_constraints, 
                                            const HypernodeID max_constraints_per_node) {
    // Read Hypergraph
    HyperedgeID num_edges = 0;
    HypernodeID num_nodes = 0;
    HyperedgeID num_removed_single_pin_hyperedges = 0;
    io::HyperedgeVector hyperedges;
    vec<HyperedgeWeight> hyperedges_weight;
    vec<HypernodeWeight> hypernodes_weight;

    io::readHypergraphFile(hg_path.string(), num_edges, num_nodes, num_removed_single_pin_hyperedges,
                         hyperedges, hyperedges_weight, hypernodes_weight);
    ALWAYS_ASSERT(hyperedges.size() == num_edges);
    ALWAYS_ASSERT(num_removed_single_pin_hyperedges == 0);

    if (num_constraints <= 0) {
        num_constraints = num_nodes * DEFAULT_CONSTRAINT_FRACTION;
    }

    vec<vec<NodeID>> adjacency;
    adjacency.resize(num_nodes);
    HypernodeID generated_constraints = get_constraints(adjacency, hyperedges, max_constraints_per_node, num_constraints);
    std::ofstream out_stream(constraints_path);
    for (NodeID i = 0; i < num_nodes; i++) {
        vec<NodeID> neighbors = adjacency[i];
        if(!neighbors.empty()) {
            std::sort(neighbors.begin(), neighbors.end());
            for(NodeID node : neighbors){
                out_stream << i << " " << node << std::endl;
            }
        }
        
    }
    out_stream.close();
    return generated_constraints;
}

HypernodeID generate_constraints_from_partitioned_hg(const fs::path hg_path, 
                                                        const fs::path part_hg_path, 
                                                        const fs::path constraints_path, 
                                                        HypernodeID num_constraints, 
                                                        const HypernodeID max_constraints_per_node) {
    HyperedgeID num_edges;
    HypernodeID num_nodes;
    std::vector<PartitionID> partitions;
    io::onlyReadHGRHeader(hg_path.string(), num_edges, num_nodes);
    io::readPartitionFile(part_hg_path.string(), num_nodes, partitions);

    if (num_constraints <= 0) {
        num_constraints = num_nodes * DEFAULT_CONSTRAINT_FRACTION;
    }

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
    num_constraints = std::min(num_constraints, constraint_count);
    print_constraints(constraints_path, constraint_list, num_constraints);
    logNumberOfDifferentNodes(constraints_path, constraint_list, num_constraints);
    return num_constraints;
}

HypernodeID generate_fixed_from_partitioned_hg(const fs::path hg_path, 
                                                        const fs::path part_hg_path, 
                                                        const fs::path fixed_path, 
                                                        HypernodeID num_constraints, 
                                                        const HypernodeID max_constraints_per_node) {
    HyperedgeID num_edges;
    HypernodeID num_nodes;
    std::vector<PartitionID> partitions;
    io::onlyReadHGRHeader(hg_path.string(), num_edges, num_nodes);
    io::readPartitionFile(part_hg_path.string(), num_nodes, partitions);

    if (num_constraints <= 0) {
        num_constraints = num_nodes * DEFAULT_CONSTRAINT_FRACTION;
    }

    std::ofstream out_stream(fixed_path);
    HypernodeID count = 0;
    HypernodeID fixed = num_nodes / 30565;
    for (HypernodeID i = 0; i < num_nodes; i++) {
        PartitionID part = -1;
        if (i%fixed == 0){ 
            part = partitions[i];
            count++;
        }
        out_stream << part << std::endl;
    }
    out_stream.close();
    LOG << "Generated fixed vertecies:" << count;
    return count;
}

int main(int argc, char* argv[]) {
    fs::path hypergraph_path;
    vec<fs::path> hypergraph_files;
    std::optional<fs::path> partitioned_hypergraph_path;
    fs::path constraint_dir;
    HypernodeID k;
    HypernodeID num_constraints;
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
        po::value<HypernodeID>(&num_constraints)->value_name("<int>")->default_value(0),
        "Number of constraints (optional)")
        ("part-hypergraph,p",
        po::value<fs::path>()->notifier([&](const fs::path& path) {
            partitioned_hypergraph_path = path;
        })->value_name("<path>"),
        "Partitioned-Hypergraph Filename or directory (optional)")
    ;
    po::variables_map cmd_vm;
    po::store(po::parse_command_line(argc, argv, options), cmd_vm);
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
            generated_constraints = generate_constraints_from_partitioned_hg(hg_file, part_hg_file.value(), constraint_file, num_constraints, max_constraints_per_node);
        } else {
            generated_constraints = generate_constraints_from_hg(hg_file, constraint_file, num_constraints, max_constraints_per_node);
        }
        LOG << "";
        LOG << "Generated " << generated_constraints << "constraints for hg:" << hg_file.filename();
    }
    LOG << "";
    return 0;
}
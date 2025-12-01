#include <boost/program_options.hpp>

#include <optional>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>


#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"

using namespace mt_kahypar;
namespace po = boost::program_options;
namespace fs = std::filesystem;

constexpr double DEFAULT_CONSTRAINT_FRACTION = 0.25;
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
                                            std::ofstream& out_stream, 
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
    for (NodeID i = 0; i < num_nodes; i++) {
        vec<NodeID> neighbors = adjacency[i];
        if(!neighbors.empty()) {
            std::sort(neighbors.begin(), neighbors.end());
            for(NodeID node : neighbors){
                out_stream << i << " " << node << std::endl;
            }
        }
        
    }
    return generated_constraints;
}

HypernodeID generate_constraints_from_partitioned_hg(const fs::path hg_path, 
                                                        const fs::path part_hg_path, 
                                                        std::ofstream& out_stream, 
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

    HypernodeID constraint_count = 0;
    std::unordered_map<HypernodeID, HypernodeID> constraints_per_node;
    for (HypernodeID node = 0; node < num_nodes; node++) {
        PartitionID node_partition = partitions[node];
        for (HypernodeID other_node = node + 1; other_node < num_nodes; other_node++) {
            if (constraint_count >= num_constraints){
                    return constraint_count;
                }
            if (constraints_per_node[node] >= max_constraints_per_node) break;
            if (constraints_per_node[other_node] >= max_constraints_per_node) continue;

            if (node_partition == partitions[other_node]) {
                constraint_count++;
                constraints_per_node[node]++;
                constraints_per_node[other_node]++;
                out_stream << node << " " << other_node << std::endl;
            }
        }
        constraints_per_node.erase(node);
    }
    return constraint_count;
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
        std::ofstream out_stream(constraint_file);
        HypernodeID generated_constraints;
        if (partitioned_hypergraph_path) {
            auto part_hg_file = findFileWithPrefix(partitioned_hypergraph_path.value(), hg_file.filename().string());
            if (!part_hg_file){
                throw fs::filesystem_error(
                    "no matching partitioned hypergraph file found in ",
                    partitioned_hypergraph_path.value(),
                    std::make_error_code(std::errc::not_a_directory)
                );
            }
            generated_constraints = generate_constraints_from_partitioned_hg(hg_file, part_hg_file.value(), out_stream, num_constraints, max_constraints_per_node);
        } else {
            generated_constraints = generate_constraints_from_hg(hg_file, out_stream, num_constraints, max_constraints_per_node);
        }
        out_stream.close();

        LOG << "";
        LOG << "Generated " << generated_constraints << "constraints for hg:" << hg_file.filename();
        LOG << "";
    }
    
    return 0;
}
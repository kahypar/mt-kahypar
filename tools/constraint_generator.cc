#include <boost/program_options.hpp>

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

HypernodeID get_constraints(vec<vec<NodeID>>& adjacency, const io::HyperedgeVector& hyperedges, const HypernodeID max_constraints_per_node, const HypernodeID num_constraints) {
    HypernodeID constraint_count = 0;
    std::unordered_map<HypernodeID, HypernodeID> constraint_per_node;
    for (io::Hyperedge edge : hyperedges) {
        for (HyperedgeID i = 0; i < edge.size(); i++) {
            HypernodeID node = edge[i];
            for (HyperedgeID j = i + 1; j < edge.size(); j++) {
                HypernodeID other_node = edge[j];
                if (constraint_count >= num_constraints){
                    return constraint_count;
                }
                if (constraint_per_node[node] >= max_constraints_per_node) break;
                if (constraint_per_node[other_node] >= max_constraints_per_node) continue;
                
                auto& list = adjacency[node];
                // if constraint already exists
                if (std::find(list.begin(), list.end(), other_node) != list.end()) continue;
                auto& other_list = adjacency[other_node];
                if (std::find(other_list.begin(), other_list.end(), node) != other_list.end()) continue;
                
                constraint_count++;
                constraint_per_node[node]++;
                constraint_per_node[other_node]++;
                adjacency[node].push_back(other_node);
            }
        }
    }
    return constraint_count;
}

void generate_constraints_for_hg(const fs::path hg_path, std::ofstream& out_stream, HypernodeID num_constraints, const HypernodeID max_constraints_per_node) {
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

    LOG << "";
    LOG << "Generated " << generated_constraints << "constraints for hg:" << hg_path.filename();
    LOG << "";
}

int main(int argc, char* argv[]) {
    fs::path hypergraph_path;
    vec<fs::path> hypergraph_files;
    fs::path constraint_dir;
    HypernodeID k;
    HypernodeID num_constraints;
    HypernodeID max_constraints_per_node;

    po::options_description options("Options");
    options.add_options()
        ("hypergraph,h",
        po::value<fs::path>(&hypergraph_path)->value_name("<path>")->required(),
        "Hypergraph Filename or hypergraph directory")
        ("constraint,c",
        po::value<fs::path>(&constraint_dir)->value_name("<path>")->required(),
        "Constraint directory")
        ("blocks,k",
        po::value<HypernodeID>(&k)->value_name("<int>")->required(),
        "Number of blocks")
        ("num-constraints,n",
        po::value<HypernodeID>(&num_constraints)->value_name("<int>")->default_value(0),
        "Number of constraints (optional)")
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

    if (fs::is_regular_file(hypergraph_path)) {
        hypergraph_files.push_back(hypergraph_path);
    } else if (fs::is_directory(hypergraph_path)) {
        for (const auto& file : fs::directory_iterator(hypergraph_path)) {
            if (fs::is_regular_file(file.path())) {
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
        generate_constraints_for_hg(hg_file, out_stream, num_constraints, max_constraints_per_node);
        out_stream.close();
    }
    
    return 0;
}
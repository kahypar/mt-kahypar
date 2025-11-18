#include <boost/program_options.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>


#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"

using namespace mt_kahypar;
namespace po = boost::program_options;

constexpr double DEFAULT_CONSTRAINT_FRACTION = 0.15;

HypernodeID constraints(std::ofstream& out_stream, const io::HyperedgeVector& hyperedges, const PartitionID& max_constraints_per_node, const HypernodeID& num_constraints) {
    HypernodeID constraint_count = 0;
    std::unordered_map<HypernodeID, HypernodeID> constraint_per_node;
    for (io::Hyperedge edge : hyperedges) {
        for (HypernodeID i = 0; i < edge.size(); i++) {
            HypernodeID node = edge[i];
            for (HyperedgeID j = i + 1; j < edge.size(); j++) {
                if (constraint_count >= num_constraints){
                    return constraint_count;
                }
                if (constraint_per_node[node] >= max_constraints_per_node) {
                    break;
                }
                HypernodeID other_node = edge[j];
                if (constraint_per_node[other_node] >= max_constraints_per_node) {
                    continue;
                }
                constraint_count++;
                constraint_per_node[node]++;
                out_stream << node << " " << other_node << std::endl;
            }
        }
    }
    return constraint_count;
}

int main(int argc, char* argv[]) {
    std::string hypergraph_file;
    std::string constraint_file;
    PartitionID max_constraints_per_node;
    HypernodeID num_constraints;

    po::options_description options("Options");
    options.add_options()
        ("hypergraph,h",
        po::value<std::string>(&hypergraph_file)->value_name("<string>")->required(),
        "Hypergraph Filename")
        ("constraint,c",
        po::value<std::string>(&constraint_file)->value_name("<string>")->required(),
        "Constraint Filename")
        ("blocks,k",
        po::value<PartitionID>(&max_constraints_per_node)->value_name("<int>")->required(),
        "Number of blocks")
        ("num-constraints,n",
        po::value<HypernodeID>(&num_constraints)->value_name("<int>")->default_value(0),
        "Number of constraints (optional)")
    ;
    po::variables_map cmd_vm;
    po::store(po::parse_command_line(argc,argv, options), cmd_vm);
    po::notify(cmd_vm);
    max_constraints_per_node--; // max constraints per node are #Blocks -1

    std::ofstream out_stream(constraint_file.c_str());
    // Read Hypergraph
    HyperedgeID num_edges = 0;
    HypernodeID num_nodes = 0;
    HyperedgeID num_removed_single_pin_hyperedges = 0;
    io::HyperedgeVector hyperedges;
    vec<HyperedgeWeight> hyperedges_weight;
    vec<HypernodeWeight> hypernodes_weight;

    io::readHypergraphFile(hypergraph_file, num_edges, num_nodes, num_removed_single_pin_hyperedges,
                         hyperedges, hyperedges_weight, hypernodes_weight);
    ALWAYS_ASSERT(hyperedges.size() == num_edges);
    ALWAYS_ASSERT(num_removed_single_pin_hyperedges == 0);

    if (num_constraints <= 0) {
        num_constraints = num_nodes * DEFAULT_CONSTRAINT_FRACTION;
    }

    LOG << "";
    LOG << "Generated " << constraints(out_stream, hyperedges, max_constraints_per_node, num_constraints) << "constraints";
    LOG << "";

    out_stream.close();
    return 0;
}
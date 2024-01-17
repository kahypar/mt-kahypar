#!/usr/bin/python3
import argparse

# Converts a graph (in metis format) into a hypergraph (in hmetis format)
# by interpreting nodes as hyperedges and edges as nodes. A hyperedge
# connects all edges of the original graph that are incident to the same node.

# Each node in the resulting hypergraph has degree 2. Further, a vertex
# partition of the hypergraph corresponds to an edge partition of the original
# graph.

parser = argparse.ArgumentParser()
parser.add_argument("input", type=str)
parser.add_argument("hypergraph", type=str)

args = parser.parse_args()

mapping = {}

def insertToMapping(u):
    if u in mapping:
        mapping[u] = mapping[u] + 1
    else:
        mapping[u] = 1
    return u

with open(args.hypergraph, "w") as out:
    with open(args.input) as f:
        header = f.readline()
        while header.startswith("%"):
            header = f.readline()
        header_vals = header.strip().split(" ")
        n_nodes = int(header_vals[0])
        n_edges = int(header_vals[1])
        dim = int(header_vals[3])
        out.write(' '.join([header_vals[1], header_vals[0], header_vals[2]]))
        i = 0
        for line in f:
            vals = line.strip().split(" ")
            i = i + 1
            x = list(filter(lambda x : int(x) > i, vals[dim:]))
            if len(x) == 0:
                continue
            out.write("\n" + str(i) + " ")
            out.write(("\n" + str(i) + " ").join(x))
        out.write("\n")
    with open(args.input) as f:
        header = f.readline()
        print(header)
        while header.startswith("%"):
            header = f.readline()
        header_vals = header.strip().split(" ")
        n_nodes = int(header_vals[0])
        n_edges = int(header_vals[1])
        dim = int(header_vals[3])
        for line in f:
            vals = line.strip().split(" ")
            out.write(' '.join([vals[0], vals[1]]))
            out.write("\n")
        
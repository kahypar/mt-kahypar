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

with open(args.input) as f, open(args.hypergraph, "w") as out:
    header = f.readline()
    while header.startswith("%"):
        header = f.readline()
    header_vals = header.strip().split(" ")
    n_edges = int(header_vals[0])
    n_nodes = int(header_vals[1])
    assert len(header_vals) <= 3
    out.write(' '.join([header_vals[0], header_vals[1], "11\n"]))
    counter = n_edges
    for line in f:
        if counter == 0:
            break
        counter = counter -1
        vals = line.strip().split(" ")
        out.write(vals[0])
        out.write(" ")
        out.write(' '.join([str(insertToMapping(u)) for u in map(int, vals[1:])]))
        out.write("\n")
    if header_vals[2] == "11":
        for line in f:
            out.write(line[:-1])
            out.write(" 1\n")
    else:
        for i in range(1, n_nodes + 1):
            if i in mapping:
                out.write(str(mapping[i]))
            else:
                out.write("0")
            out.write(" 1\n")
# File Formats

The following is an overview of the input and output file formats which are used by Mt-KaHyPar.

**Important note:** For historical reasons, the hMetis and Metis input formats use 1-based indexing.
However, Mt-KaHyPar converts the indices to 0-based indexing when reading the files.
Any results obtained from the binary or one of the library interfaces will also use 0-based indices (i.e., shifted by -1 in comparison to the input file).

### hMetis Format for Input Hypergraphs

Per default, Mt-KaHyPar assumes that the input hypergraph is provided in the hMetis format: [Example](tests/instances/unweighted_hypergraph.hgr)

The general format looks as follows:

```
% comment
num_hypergedges num_hypernodes [weight_type]
[hyperedge_weight] pin_1 pin_2 ... pin_i
...
[hyperedge_weight] pin_1 pin_2 ... pin_j
[hypernode_weight_1]
...
[hypernode_weight_n]
```

Any line that starts with ‘%’ is a comment line and is skipped.
The first line is a header containing two or three numbers describing the total number of hyperedges, the total number of hypernodes and the types weights used by the hypergraph
(00/omitted = unweighted, 10 = node weight, 01 = edge weight, 11 = node and edge weights).

Afterwards, there is one line for each hyperedge which contains a list of the pins (hypernode IDs) of the hyperedge.
If hyperedge weights are used, there is an additional entry at the start of the line which is the weight of the hyperedge.

If no hypernode weights are used, this is the end of the file.
Otherwise, there is one line for each hypernode containing a single entry, which is the weight of the hypernode.

#!/usr/bin/python3
import random
input = "/home/konstantin/hypergraph_partitioner/ibm01.weighted.metis"
output = "/home/konstantin/hypergraph_partitioner/benchmark_set/test.metis"
before_dims=1
new_dims=2
with open(input) as ip:
    with open(output, "w") as op:
        data = ip.readlines()
        op.write(data[0])
        for i in range(1, len(data)):
            splitted = data[i].split(" ")
            res = splitted[0:before_dims]
            for j in range(0,new_dims):
                res = res + [str(random.randrange(0,len(data)))]

            op.write(' '.join(res + splitted[before_dims:]))
             

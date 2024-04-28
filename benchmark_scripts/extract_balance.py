#!/usr/bin/python3
input = "/home/konstantin/hypergraph_partitioner/04_27_metis_compare.txt"
with open(input) as ip:
    counter = 0
    total = 0
    this_counter = 0
    for line in ip:
        if "#0" in line or "#1" in line or "#2" in line:
            s = float(line.strip().split(":")[1].strip().split(" ")[0])
            if s > 1.06:
                this_counter = 1
        if "#2" in line:
            counter = counter + this_counter
            total = total + 1
            this_counter = 0
    print(str(counter / total))
print(counter)
print(total)
#!/usr/bin/python3
import subprocess
import ntpath
import argparse
import time
import re
import math
import os
import os.path
import shutil
from threading import Timer
import signal

###################################
# SETUP ENV
###################################
algorithm = "hMetis-R"
hmetis = os.environ.get("HMETIS")
assert (hmetis != None), "check env.sh"
###################################

parser = argparse.ArgumentParser()
parser.add_argument("graph", type=str)
parser.add_argument("k", type=int)
parser.add_argument("epsilon", type=float)
parser.add_argument("seed", type=int)
parser.add_argument("objective", type=str)
parser.add_argument("timelimit", type=int)
parser.add_argument("--partition_folder", type=str, default = "")
parser.add_argument("--config", type=str, default = "")
parser.add_argument("--name", type=str, default = "")

args = parser.parse_args()

if args.name != "":
  algorithm = args.name

# Read Number of Nodes
hg = open(args.graph, 'r')
for line in hg:
     # ignore comment lines
    if line.startswith('%'):
        continue
    hg_param = line.split()
    numNodes= hg_param[1]
    break;
hg.close()

if args.objective == "cut":
  objective = "cut"
elif args.objective == "km1":
  objective = "soed"

graph_file = args.graph
if args.partition_folder != "":
   graph_file = args.partition_folder + "/" + ntpath.basename(args.graph) + ".part" + str(args.k) + ".epsilon" + str(args.epsilon) + ".seed" + str(args.seed)
   shutil.copyfile(args.graph, graph_file)

mode = "rb"
#We use hMetis-RB as initial partitioner. If called to partition a graph into k parts
#with an UBfactor of b, the maximal allowed partition size will be 0.5+(b/100)^(log2(k)) n.
#In order to provide a balanced initial partitioning, we determine the UBfactor such that
#the maximal allowed partiton size corresponds to our upper bound i.e.
#(1+epsilon) * ceil(total_weight / k).
exp = 1.0 / math.log(args.k,2)
ufactor = 50.0 * (2 * math.pow((1 + args.epsilon), exp)
          * math.pow(math.ceil(float(numNodes)/args.k) / float(numNodes), exp) - 1)
if ufactor < 0.1:
  ufactor = 0.1

hmetis_command = [hmetis,
                  str(graph_file),
                  str(args.k),
                  '-ptype=' + mode,
                  '-dbglvl=34',
                  '-otype=' + objective,
                  '-ufactor='+str(ufactor),
                  '-seed='+str(args.seed)]
if objective == "soed":
  hmetis_command.extend(["-reconst"])
hmetis_proc = subprocess.Popen(hmetis_command, stdout=subprocess.PIPE, universal_newlines=True, preexec_fn=os.setsid)

def kill_proc():
	os.killpg(os.getpgid(hmetis_proc.pid), signal.SIGTERM)

t = Timer(args.timelimit, kill_proc)
t.start()
out, err = hmetis_proc.communicate()
t.cancel()
end = time.time()

part_sizes = []
total_time = 2147483647
cut = 2147483647
km1 = 2147483647
soed = 2147483647
imbalance = 1.0
timeout = "no"
failed = "no"

if hmetis_proc.returncode == 0:
  # Extract metrics out of hMetis output
  for line in out.split('\n'):
    s = str(line).strip()
    if ("Hyperedge Cut" in s):
      cut = int(s.split()[2].split('.')[0])
    if ("Sum of External" in s):
      soed = int(s.split()[4].split('.')[0])
      km1 = soed - cut
    if ("Multilevel" in s):
      total_time = float(s.split()[1])
    if ("[" in s):
      split = s.split(']')
      for token in split[:-1]:
        t = re.search('\(.*?\)', token)
        part_size = float(t.group(0)[1:-1])
        part_sizes.append(part_size)
elif hmetis_proc.returncode == -signal.SIGTERM:
  timeout = "yes"
else:
  failed = "yes"

if len(part_sizes) > 0:
  # compute imbalance
  max_part_size = max(part_sizes)
  total_weight = sum(part_sizes)
  imbalance = float(max_part_size) / math.ceil(float(total_weight)/args.k) - 1.0

# CSV format: algorithm,graph,timeout,seed,k,epsilon,num_threads,imbalance,totalPartitionTime,objective,km1,cut,failed
print(algorithm,
      ntpath.basename(args.graph),
      timeout,
      args.seed,
      args.k,
      args.epsilon,
      1,
      imbalance,
      total_time,
      args.objective,
      km1,
      cut,
      failed,
      sep=",")

if args.partition_folder != "":
   os.remove(graph_file)
   src_partition_file = graph_file + ".part." + str(args.k)
   dst_partition_file = args.partition_folder + "/" + ntpath.basename(args.graph) + ".part" + str(args.k) + ".epsilon" + str(args.epsilon) + ".seed" + str(args.seed) + ".partition"
   if os.path.exists(src_partition_file):
    shutil.move(src_partition_file, dst_partition_file)

#!/usr/bin/python3
import subprocess
import ntpath
import argparse
import time
import re
import math
import os
import os.path
from threading import Timer
import signal

###################################
# SETUP ENV
###################################
algorithm = "Metis-R"
metis = os.environ.get("METIS")
assert (metis != None), "check env.sh"
###################################

parser = argparse.ArgumentParser()
parser.add_argument("graph", type=str)
parser.add_argument("k", type=int)
parser.add_argument("epsilon", type=float)
parser.add_argument("seed", type=int)
parser.add_argument("objective", type=str)
parser.add_argument("timelimit", type=int)
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
    numNodes= hg_param[0]
    break;
hg.close()

#We use hMetis-RB as initial partitioner. If called to partition a graph into k parts
#with an UBfactor of b, the maximal allowed partition size will be 0.5+(b/100)^(log2(k)) n.
#In order to provide a balanced initial partitioning, we determine the UBfactor such that
#the maximal allowed partiton size corresponds to our upper bound i.e.
#(1+epsilon) * ceil(total_weight / k).
exp = 1.0 / math.log(args.k,2)
ufactor = (math.pow((1 + args.epsilon), exp) - 1.0) * 1000.0

metis_command = [metis,
                 str(args.graph),
                 str(args.k),
                 '-ptype=rb',
                 '-dbglvl=1',
                 '-objtype=cut',
                 '-ufactor='+str(ufactor),
                 '-nooutput',
                 '-seed='+str(args.seed)]

metis_proc = subprocess.Popen(metis_command, stdout=subprocess.PIPE, universal_newlines=True, preexec_fn=os.setsid)

def kill_proc():
	os.killpg(os.getpgid(metis_proc.pid), signal.SIGTERM)

t = Timer(args.timelimit, kill_proc)
t.start()
out, err = metis_proc.communicate()
t.cancel()
end = time.time()

total_time = 2147483647
cut = 2147483647
km1 = 2147483647
soed = 2147483647
imbalance = 1.0
timeout = "no"
failed = "no"

if metis_proc.returncode == 0:
  # Extract metrics out of hMetis output
  for line in out.split('\n'):
    s = str(line).strip()
    if ("Edgecut:" in s):
      cut = int(s.split()[2].split(',')[0])
      km1 = cut
    if ("Partitioning:" in s):
      total_time = float(s.split()[1])
    if ("actual:" in s):
      max_part_size = float(s.split()[3].split(',')[0])
      imbalance = max_part_size / math.ceil(float(numNodes)/args.k) - 1.0
elif metis_proc.returncode == -signal.SIGTERM:
  timeout = "yes"
else:
  failed = "yes"

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

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
algorithm = "Zoltan"
zoltan = os.environ.get("ZOLTAN")
experiment_dir = os.path.join(os.getcwd(), 'zoltan_driver_temps')
assert (zoltan != None), "check env.sh"
###################################

parser = argparse.ArgumentParser()
# make sure it's in Zoltan's format and ends on .zoltan.hg
# same as PaToH's format but every node weight on its own line.
parser.add_argument("graph", type=str)
parser.add_argument("threads", type=int)
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

graph_absolute_path = os.path.abspath(args.graph)
graph_path_without_extension = "".join(graph_absolute_path.rsplit(".hg", 1))
graph_name = ntpath.basename(args.graph.replace(".zoltan.hg", ""))

# Zoltan takes its parameters via a file 'zdrive.inp' in the current
# working directory --> create a uniqe wd for each call
wd = experiment_dir + "/" + ntpath.basename(graph_path_without_extension) + "/thread" + str(args.threads) + "/k" + str(args.k) + "/seed" + str(args.seed)
os.makedirs(wd, exist_ok=True)
os.chdir(wd)

# write weird Zoltan input file
with open('zdrive.inp', 'w') as f:
	print("Decomposition Method = hypergraph", file=f)
	print("Zoltan Parameters = HYPERGRAPH_PACKAGE=phg", file=f)
	print("Zoltan Parameters = DEBUG_LEVEL=0", file=f)
	print("Zoltan Parameters = LB_APPROACH=PARTITION", file=f)
	if args.objective == "cut":
		print("Zoltan Parameters = PHG_CUT_OBJECTIVE=HYPEREDGES", file=f)
	elif args.objective == "km1":
		print("Zoltan Parameters = PHG_CUT_OBJECTIVE=CONNECTIVITY", file=f)
	else:
		raise RuntimeError("Unknown objective: " + args.objective)
	print("Zoltan Parameters = NUM_GLOBAL_PARTITIONS=", args.k, file=f)
	print("Zoltan Parameters = IMBALANCE_TOL=", 1 + args.epsilon, file=f)
	print("Zoltan Parameters = PHG_EDGE_SIZE_THRESHOLD=1.0", file=f)
	print("Zoltan Parameters = SEED=", args.seed, file=f)
	print("File Type = hypergraph", file=f)
	print("Compression = uncompressed", file=f)
	print("text output = 0", file=f)
	print("File Name = ", graph_path_without_extension, file=f)

# Run Zoltan
zoltan_proc = subprocess.Popen(["mpirun -N " + str(args.threads) + " " + zoltan],
                               stdout=subprocess.PIPE, universal_newlines=True, shell=True, preexec_fn=os.setsid)

def kill_proc():
  os.killpg(os.getpgid(zoltan_proc.pid), signal.SIGTERM)

t = Timer(args.timelimit, kill_proc)
t.start()
out, err = zoltan_proc.communicate()
t.cancel()
end = time.time()

total_time = 2147483647
cut = 2147483647
km1 = 2147483647
imbalance = 1.0
timeout = "no"
failed = "no"

if zoltan_proc.returncode == 0:
  i = 0
  j = 0
  k = 0
  for line in out.split('\n'):
    s = str(line).strip()
    if "CUTN" in s:
      if (i == 0):
        i=1
        continue
      cut = int(float(s.split(':')[1]))
    if "CUTL " in s:
      if (j == 0):
        j=1
        continue
      km1 = int(float(s.split(':')[1]))
    if "Zoltan_LB_Partition" in s:
      total_time = float(s.split('=')[1])
    if ("Zoltan_LB_Eval_HG  Number of objects : " in s):
        if (k == 0):
          k=1
          continue
        res = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", s)
        imbalance = float(res[3]) - 1.0

elif zoltan_proc.returncode == -signal.SIGTERM:
  timeout = "yes"
else:
  failed = "yes"

# CSV format: algorithm,graph,timeout,seed,k,epsilon,num_threads,imbalance,totalPartitionTime,objective,km1,cut,failed
print(algorithm,
      graph_name,
      timeout,
      args.seed,
      args.k,
      args.epsilon,
      args.threads,
      imbalance,
      total_time,
      args.objective,
      km1,
      cut,
      failed,
      sep=",")

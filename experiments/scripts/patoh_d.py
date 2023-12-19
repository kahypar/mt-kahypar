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
algorithm = "PaToH-D"
preset = "D"
patoh = os.environ.get("PATOH")
assert (patoh != None), "check env.sh"
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

if args.objective == "cut":
  objective = "U"
elif args.objective == "km1":
  objective = "O"

# Read hypergraph weight
with open(str(args.graph)) as hypergraph:
  header_parsed = False
  is_weighted = False
  total_weight = 0
  for line in hypergraph:
    if header_parsed or line.startswith('%'):
      continue
    else:
      header_parsed = True
      hg_params = line.split()
      # patoh file format uses 0,1,2,3 for weight types
      is_weighted = len(hg_params) >= 5 and hg_params[4] in ['1', '3']
      if not is_weighted:
        total_weight = int(hg_params[1])
        break

  if is_weighted:
    # read weight from last line
    last_line = line
    for weight in line.split():
      w = int(weight.strip())
      total_weight += w

# Run PaToH-S
patoh_proc = subprocess.Popen([patoh,
                               args.graph,
                               str(args.k),
                               'FI='+str(args.epsilon), # imbalance ratio
                               'PQ=' + preset, # preset
                               'SD=' + str(args.seed), # seed
                               'OD=1', # non verbose output
                               'UM=' + objective, # connectivity metric
                               'WI=0', # dont write the partitioning info to disk
                               'BO=C', # balance on cell
                               'A0=100', #  MemMul_CellNet
                               'A1=100', #  MemMul_Pins
                               'A2=100', #  MemMul_General
                              ], stdout=subprocess.PIPE, universal_newlines=True, preexec_fn=os.setsid)

def kill_proc():
	os.killpg(os.getpgid(patoh_proc.pid), signal.SIGTERM)

t = Timer(args.timelimit, kill_proc)
t.start()
out, err = patoh_proc.communicate()
t.cancel()
end = time.time()

total_time = 2147483647
cut = 2147483647
km1 = 2147483647
imbalance = 1.0
timeout = "no"
failed = "no"

if patoh_proc.returncode == 0:
  # Extract metrics out of PaToH output
  for line in out.split('\n'):
    s = str(line).strip()
    if ("Cells" in s):
      t = re.compile('Cells : \s*([^\s]*)')
      numHNs = int(str(t.findall(s)[0]))
    if ("'Con - 1' Cost" in s):
      km1 = int(s.split("'Con - 1' Cost:")[1])
    if ("Cut Cost" in s):
      cut = int(s.split("Cut Cost:")[1])
    if ("Part Weights" in s):
      t = re.compile('Min=\s*([^\s]*)')
      min_part = float(t.findall(s)[0])
      t = re.compile('Max=\s*([^\s]*)')
      max_part = float(t.findall(s)[0])
      imbalance = float(max_part) / math.ceil(float(total_weight) / args.k) - 1.0
    if ("Total   " in s):
      t = re.compile('Total\s*:\s*([^\s]*)')
      total_time = float(t.findall(s)[0])
elif patoh_proc.returncode == -signal.SIGTERM:
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

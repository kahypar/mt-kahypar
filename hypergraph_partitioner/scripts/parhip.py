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
algorithm = "ParHIP"
parhip = os.environ.get("PARHIP")
assert (parhip != None), "check env.sh"
###################################

parser = argparse.ArgumentParser()
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

# Run ParHIP
parhip_proc = subprocess.Popen(["mpirun -N " +str(args.threads) + " " +
                                parhip + " " +
                                args.graph + " "
                                "--k=" + str(args.k) + " "
                                "--seed=" + str(args.seed) + " "
                                "--imbalance=" + str(int(args.epsilon * 100.0)) + " "
                                "--preconfiguration=fastsocial"],
                                stdout=subprocess.PIPE, universal_newlines=True, shell=True, preexec_fn=os.setsid)

def kill_proc():
	os.killpg(os.getpgid(parhip_proc.pid), signal.SIGTERM)

t = Timer(args.timelimit, kill_proc)
t.start()
out, err = parhip_proc.communicate()
t.cancel()
end = time.time()

total_time = 2147483647
cut = 2147483647
km1 = 2147483647
imbalance = 1.0
timeout = "no"
failed = "no"

if parhip_proc.returncode == 0:
  # Extract metrics out of ParHIP output
  for line in out.split('\n'):
    s = str(line).strip()
    if "final edge cut" in s:
      cut = int(s.split('final edge cut')[1])
      km1 = int(s.split('final edge cut')[1])
    if "final balance" in s:
      imbalance = float(s.split('final balance')[1]) - 1.0
    if "total partitioning time elapsed" in s:
      total_time = float(s.split('total partitioning time elapsed')[1])
elif parhip_proc.returncode == -signal.SIGTERM:
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

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
algorithm = "ParMetis"
parmetis = os.environ.get("PARMETIS")
assert (parmetis != None), "check env.sh"
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

# Run ParMetis
parmetis_proc = subprocess.Popen(["mpirun -N " +str(args.threads) + " " +
                                parmetis + " " +
                                args.graph + " "
                                "1 " + str(args.k) + " 0 0 3 " + str(args.seed) + " " + str(args.epsilon)],
                                stdout=subprocess.PIPE, universal_newlines=True, shell=True)

def kill_proc():
	os.killpg(os.getpgid(parmetis_proc.pid), signal.SIGTERM)

t = Timer(args.timelimit, kill_proc)
t.start()
out, err = parmetis_proc.communicate()
t.cancel()
end = time.time()

total_time = 2147483647
cut = 2147483647
km1 = 2147483647
imbalance = 1.0
timeout = "no"
failed = "no"

if parmetis_proc.returncode == 0:
  # Extract metrics out of hMetis output
  for line in out.split('\n'):
    s = str(line).strip()
    if ("Cut:" in s):
      cut = int(s.split('Cut:')[1].split()[0])
      km1 = cut
    if ("Total:" in s):
      total_time = float(s.split('Sum:')[1].split()[0].split(',')[0])
    if ("Balance:" in s):
      imbalance = float(s.split('Balance:')[1].split()[0]) - 1.0
elif parmetis_proc.returncode == -signal.SIGTERM:
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
      args.threads,
      imbalance,
      total_time,
      args.objective,
      km1,
      cut,
      failed,
      sep=",")
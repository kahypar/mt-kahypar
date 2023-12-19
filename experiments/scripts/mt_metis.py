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
algorithm = "Mt-Metis"
mt_metis = os.environ.get("MT_METIS")
assert (mt_metis != None), "check env.sh"
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

# Run Mt-KaHIP
mt_metis_proc = subprocess.Popen([mt_metis,
                                  "--threads=" + str(args.threads),
                                  "--balance=" + str(args.epsilon + 1.0),
                                  "--seed=" + str(args.seed),
                                  "--times",
                                  "--partstats",
                                  "--rtype=hs",
                                  args.graph,
                                  str(args.k)],
                                 stdout=subprocess.PIPE, universal_newlines=True, preexec_fn=os.setsid)

def kill_proc():
	os.killpg(os.getpgid(mt_metis_proc.pid), signal.SIGTERM)

t = Timer(args.timelimit, kill_proc)
t.start()
out, err = mt_metis_proc.communicate()
t.cancel()
end = time.time()

total_time = 2147483647
cut = 2147483647
km1 = 2147483647
imbalance = 1.0
timeout = "no"
failed = "no"

if mt_metis_proc.returncode == 0:
  for line in out.split('\n'):
    s = str(line).strip()
    if "Best Objective" in s:
      cut = int(s.split('Best Objective:')[1])
      km1 = int(s.split('Best Objective:')[1])
    if "Total Time" in s:
      total_time = float(s.split('Total Time:')[1][:-1])
    if "constraint #0" in s:
      imbalance = float(s.split("constraint #0:")[1].split(' ')[2]) - 1.0

elif mt_metis_proc.returncode == -signal.SIGTERM:
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

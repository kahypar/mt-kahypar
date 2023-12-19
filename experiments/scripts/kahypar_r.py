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
import shutil

###################################
# SETUP ENV
###################################
algorithm = "KaHyPar-R"
kahypar_r = os.environ.get("KAHYPAR")
kahypar_r_config = os.environ.get("KAHYPAR_R_CONFIG")
assert (kahypar_r != None and kahypar_r_config != None), "check env.sh"
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

if args.config != "":
  kahypar_r_config = args.config
if args.name != "":
  algorithm = args.name

# Run KaHyPar-R
kahypar_r_proc = subprocess.Popen([kahypar_r,
                                   "-h" + args.graph,
                                   "-k" + str(args.k),
                                   "-e" + str(args.epsilon),
                                   "--seed=" + str(args.seed),
                                   "-o" + str(args.objective),
                                   "-mrecursive",
                                   "-p" + kahypar_r_config,
                                   "--sp-process=true"],
                                  stdout=subprocess.PIPE, universal_newlines=True, preexec_fn=os.setsid)

def kill_proc():
	os.killpg(os.getpgid(kahypar_k_proc.pid), signal.SIGTERM)

t = Timer(args.timelimit, kill_proc)
t.start()
out, err = kahypar_r_proc.communicate()
t.cancel()
end = time.time()

total_time = 2147483647
cut = 2147483647
km1 = 2147483647
imbalance = 1.0
timeout = "no"
failed = "no"

if kahypar_r_proc.returncode == 0:
  # Extract metrics out of KaHyPar-CA output
  for line in out.split('\n'):
    s = str(line).strip()
    if "RESULT" in s:
      km1 = int(s.split(" km1=")[1].split(" ")[0])
      cut = int(s.split(" cut=")[1].split(" ")[0])
      total_time = float(s.split(" totalPartitionTime=")[1].split(" ")[0])
      imbalance = float(s.split(" imbalance=")[1].split(" ")[0])
elif kahypar_r_proc.returncode == -signal.SIGTERM:
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

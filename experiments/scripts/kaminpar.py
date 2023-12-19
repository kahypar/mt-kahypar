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
algorithm = "KaMinPar"
kaminpar = os.environ.get("KAMINPAR")
assert (kaminpar != None), "check env.sh"
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

kaminpar_call = [kaminpar,
                 "-G" + args.graph,
                 "-k" + str(args.k),
                 "--threads="+str(args.threads),
                 "--epsilon="+str(args.epsilon),
                 "--seed="+str(args.seed),
                 "--experiment"]
if args.k >= 1024:
  kaminpar_call.append("--fast-ip")

# Run KaFFPa
kaminpar_proc = subprocess.Popen(kaminpar_call,
                                stdout=subprocess.PIPE, universal_newlines=True, preexec_fn=os.setsid)

def kill_proc():
	os.killpg(os.getpgid(kaminpar_proc.pid), signal.SIGTERM)

t = Timer(args.timelimit, kill_proc)
t.start()
out, err = kaminpar_proc.communicate()
t.cancel()
end = time.time()

total_time = 2147483647
cut = 2147483647
km1 = 2147483647
imbalance = 1.0
timeout = "no"
failed = "no"

if kaminpar_proc.returncode == 0:
  for line in out.split('\n'):
    s = str(line).strip()
    if "RESULT" in s:
      cut = int(s.split(' cut=')[1].split(" ")[0])
      km1 = cut
      imbalance = float(s.split(' imbalance=')[1].split(" ")[0])
    if "TIME" in s:
      total_time = float(s.split(' partitioning=')[1].split(" ")[0])
elif kaminpar_proc.returncode == -signal.SIGTERM:
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

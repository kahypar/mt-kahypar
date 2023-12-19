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
algorithm = "Scotch"
scotch = os.environ.get("SCOTCH")
assert (scotch != None), "check env.sh"
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


scotch_command = [scotch,
                  str(args.k),
                  str(args.graph),
                  '-b'+str(args.epsilon),
                  '-crq',
                  '-vmts']


scotch_proc = subprocess.Popen(scotch_command, stdout=subprocess.PIPE, universal_newlines=True, preexec_fn=os.setsid)

def kill_proc():
  os.killpg(os.getpgid(scotch_proc.pid), signal.SIGTERM)

t = Timer(args.timelimit, kill_proc)
t.start()
out, err = scotch_proc.communicate()
t.cancel()
end = time.time()

total_time = 2147483647
cut = 2147483647
km1 = 2147483647
soed = 2147483647
imbalance = 1.0
timeout = "no"
failed = "no"

if scotch_proc.returncode == 0:
  # Extract metrics out of hMetis output
  for line in out.split('\n'):
    s = str(line).strip()
    if ("CommCutSz" in s):
      cut = int(s.split()[2].split('(')[1].split(')')[0])
      km1 = cut
    if ("Mapping" in s):
      total_time = float(s.split()[2])
    if ("Target" in s):
      imbalance = float(s.split()[6].split('=')[1]) - 1.0
elif scotch_proc.returncode == -signal.SIGTERM:
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

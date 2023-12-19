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
algorithm = "Hype"
hype = os.environ.get("HYPE")
assert (hype != None), "check env.sh"
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

# Run HYPE
hype_proc = subprocess.Popen([hype,
                              args.graph,
                              str(args.k)],
                              stdout=subprocess.PIPE, universal_newlines=True, preexec_fn=os.setsid)

def kill_proc():
	os.killpg(os.getpgid(hype_proc.pid), signal.SIGTERM)

t = Timer(args.timelimit, kill_proc)
t.start()
out, err = hype_proc.communicate()
t.cancel()
end = time.time()

total_time = 2147483647
cut = 2147483647
km1 = 2147483647
imbalance = 1.0
timeout = "no"
failed = "no"

if hype_proc.returncode == 0:
  # Extract metrics out of HYPE output
  for line in out.split('\n'):
    s = str(line).strip()
    if ("Hyperedges cut:" in s):
        t = re.compile('Hyperedges cut: \s*([^\s]*)')
        cut = int(t.findall(s)[0])
    if ("K-1:" in s):
        t = re.compile('K-1: \s*([^\s]*)')
        km1 = int(t.findall(s)[0])
    if ("node balancing:" in s):
        t = re.compile('node balancing: \s*([^\s]*)')
        imbalance = float(t.findall(s)[0])
    if ("partition time:" in s):
        t = re.compile('partition time: \s*([^\s]*)')
        total_time = float(t.findall(s)[0])/1000

elif hype_proc.returncode == -signal.SIGTERM:
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

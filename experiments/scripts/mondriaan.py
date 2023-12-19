#!/usr/bin/python3
import subprocess
import ntpath
import argparse
import time
import re
import math
import os
import os.path
import glob
from threading import Timer
import signal

###################################
# SETUP ENV
###################################
algorithm = "Mondriaan"
mondriaan = os.environ.get("MONDRIAAN")
evaluator = os.environ.get("MONDRIAAN_EVALUATOR")
assert (mondriaan != None and evaluator != None), "check env.sh"
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
  objective = "cutnet"
elif args.objective == "km1":
  objective = "lambda1"

# Run Mondriaan
start = time.time()
mondriaan_proc = subprocess.Popen([mondriaan,
                                   str(args.graph),
                                   str(args.k),
                                   str(args.epsilon),
                                   '-Metric=' + objective,
	                                 '-SplitStrategy=onedimcol',
                                   '-Seed='+str(args.seed)],
                                   stdout=subprocess.PIPE, universal_newlines=True, preexec_fn=os.setsid)

def kill_proc():
	os.killpg(os.getpgid(mondriaan_proc.pid), signal.SIGTERM)

t = Timer(args.timelimit, kill_proc)
t.start()
out, err = mondriaan_proc.communicate()
t.cancel()
end = time.time()

total_time = end - start
cut = 2147483647
km1 = 2147483647
imbalance = 1.0
timeout = "no"
failed = "no"

hgr_file = re.sub('.mondriaan.mtx$', '', args.graph)
mondriaan_output_file = args.graph+'-v'+str(args.k)+'-s'+str(args.seed)
if mondriaan_proc.returncode == 0:
  out, err = subprocess.Popen([evaluator,
                               hgr_file,
                               mondriaan_output_file],
                              stdout=subprocess.PIPE, universal_newlines=True).communicate()

  # Extract metrics out of Mondriaan Evaluator output
  for line in out.split('\n'):
    s = str(line).strip()
    if "cut" in s and "RESULT" not in s:
      cut = int(s.split('=')[1])
    if "km1" in s and "RESULT" not in s:
      km1 = int(s.split('=')[1])
    if "imbalance" in s and "RESULT" not in s:
      imbalance = float(s.split('=')[1])
elif mondriaan_proc.returncode == -signal.SIGTERM:
  timeout = "yes"
  total_time = 2147483647
else:
  failed = "yes"
  total_time = 2147483647

# CSV format: algorithm,graph,timeout,seed,k,epsilon,num_threads,imbalance,totalPartitionTime,objective,km1,cut,failed
print(algorithm,
      ntpath.basename(hgr_file),
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

def remove_if_exists(file):
  if os.path.exists(file):
    os.remove(file)

remove_if_exists(args.graph + "-local" + str(args.k) + "-s" + str(args.seed) )
remove_if_exists(args.graph + "-P" + str(args.k) + "-s" + str(args.seed) )
remove_if_exists(args.graph + "-reor-P" + str(args.k) + "-s" + str(args.seed) )
remove_if_exists(args.graph + "-v" + str(args.k) + "-s" + str(args.seed) )
remove_if_exists(args.graph + "-u" + str(args.k) + "-s" + str(args.seed) )
remove_if_exists(args.graph + "-I" + str(args.k) + "-s" + str(args.seed) )

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
import sys

parser = argparse.ArgumentParser()
parser.add_argument("partitioner", type=str)
parser.add_argument("graph", type=str)
parser.add_argument("k", type=int)
parser.add_argument("epsilon", type=float)
parser.add_argument("seed", type=int)
parser.add_argument("objective", type=str)
parser.add_argument("timelimit", type=int)

args = parser.parse_args()

seed = int(args.seed)
timelimit = float(args.timelimit)
found = False

best_objective = 2147483647
while timelimit > 0:
  partitioner_command = [args.partitioner,
                        str(args.graph),
                        str(args.k),
                        str(args.epsilon),
                        str(seed),
                        str(args.objective),
                        str(int(timelimit + 1))]
  start = time.time()
  partitioner_command = subprocess.Popen(partitioner_command, stdout=subprocess.PIPE, universal_newlines=True, preexec_fn=os.setsid)

  def kill_proc():
    os.killpg(os.getpgid(partitioner_command.pid), signal.SIGTERM)

  t = Timer(timelimit, kill_proc)
  t.start()
  out, err = partitioner_command.communicate()
  t.cancel()
  end = time.time()

  if partitioner_command.returncode == 0:
    # CSV format: algorithm,graph,timeout,seed,k,epsilon,num_threads,imbalance,
    #             totalPartitionTime,objective,km1,cut,failed
    values = out.split(',')
    partition_time = end - start # total partition time
    timelimit -= partition_time
    timeout = values[2]
    if timeout == "yes":
      timelimit = 0
    else:
      seed = seed + 10
      metric = 2147483647
      if args.objective == "km1":
        metric = int(values[10]) # connectivity metric
      elif args.objective == "cut":
        metric = int(values[11]) # cut metric

      if metric < best_objective:
        found = True
        best_objective = metric
        values[3] = str(args.seed)
        values[8] = str(float(args.timelimit) - timelimit)
        print(','.join(values).split('\n')[0])
        sys.stdout.flush()
  else:
    timelimit = 0

if not found:
  print("Algo",
      ntpath.basename(args.graph),
      "yes",
      args.seed,
      args.k,
      args.epsilon,
      1,
      1.0,
      2147483647,
      args.objective,
      2147483647,
      2147483647,
      "no",
      sep=",")



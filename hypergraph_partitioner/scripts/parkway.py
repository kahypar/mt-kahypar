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
algorithm = "Parkway"
parkway = os.environ.get("PARKWAY")
parkway_config = os.environ.get("PARKWAY_CONFIG")
hgr_to_parkway_converter = os.environ.get("HGR_TO_PARKWAY_CONVERTER")
evaluator = os.environ.get("KAHYPAR_VERIFY_PARTITION")
experiment_dir = os.path.join(os.getcwd(), 'parkway_hg_temps')
assert (parkway != None), "check env.sh"
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

if args.config != "":
  parkway_config = args.config
if args.name != "":
  algorithm = args.name

# Convert hMetis hypergraph to Parkway format
wd = experiment_dir + "/" + ntpath.basename(args.graph) + "/thread" + str(args.threads) + "/k" + str(args.k) + "/seed" + str(args.seed)
os.makedirs(wd, exist_ok=True)
parkway_file = wd + "/" + ntpath.basename(args.graph) + ".bin"
conversion_proc = subprocess.Popen([hgr_to_parkway_converter,
                                    "-h" + args.graph,
                                    "-p" + str(args.threads),
                                    "-o" + parkway_file],
                                    stdout=subprocess.PIPE, universal_newlines=True).communicate()
# Run Parkway
parkway_proc = subprocess.Popen(["mpirun -N " +str(args.threads) + " " +
                                 parkway + " " +
                                 "-p" + str(args.k) + " " +
                                 "-c" + str(args.epsilon) + " " +
                                 "--serial-partitioning.number-of-runs=1 " +
                                 "--sprng-seed=" + str(args.seed) + " " +
                                 "-o" + parkway_config + " " +
                                 "--hypergraph=" + parkway_file + " " +
                                 "--write-partitions-to-file"],
                                 stdout=subprocess.PIPE, universal_newlines=True, shell=True, preexec_fn=os.setsid)

def kill_proc():
	os.killpg(os.getpgid(parkway_proc.pid), signal.SIGTERM)

t = Timer(args.timelimit, kill_proc)
t.start()
out, err = parkway_proc.communicate()
t.cancel()
end = time.time()

total_time = 2147483647
cut = 2147483647
km1 = 2147483647
imbalance = 1.0
timeout = "no"
failed = "no"

if parkway_proc.returncode == 0:
  # Search for partition time in parkway output
  for line in out.split('\n'):
    s = str(line).strip()
    if "TOTAL TIME" in s:
      total_time = float(s.split('=')[1])

  # Evaluate Partition
  parkway_partition_file = parkway_file + ".part." + str(args.k) + "." + str(args.seed)
  evaluator_out, evaluator_err = subprocess.Popen([evaluator,
                                                  args.graph,
                                                  parkway_partition_file],
                                                  stdout=subprocess.PIPE, universal_newlines=True).communicate()

  # Extract metrics out of evaluator
  for line in evaluator_out.split('\n'):
    s = str(line).strip()
    if "cut" in s:
      cut = int(s.split('=')[1])
    if "km1" in s:
      km1 = int(s.split('=')[1])
    if "imbalance" in s:
      imbalance = float(s.split('=')[1])

elif parkway_proc.returncode == -signal.SIGTERM:
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

shutil.rmtree(wd)

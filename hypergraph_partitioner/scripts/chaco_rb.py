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
algorithm = "Chaco-R"
chaco = os.environ.get("CHACO")
assert (chaco != None), "check env.sh"
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

# Read Number of Nodes
hg = open(args.graph, 'r')
for line in hg:
     # ignore comment lines
    if line.startswith('%'):
        continue
    hg_param = line.split()
    numNodes= hg_param[0]
    break;
hg.close()

exp = 1.0 / math.log(args.k,2)
adaptive_epsilon = (math.pow((1 + args.epsilon), exp) - 1.0)
user_params = ["RANDOM_SEED = " + str(args.seed),
               "KL_IMBALANCE = " + str(adaptive_epsilon)]
with open('User_Params', 'w') as f:
  f.write('\n'.join(user_params))
chaco_command = [chaco]

chaco_proc = subprocess.Popen(chaco_command,
                              stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              universal_newlines=True, preexec_fn=os.setsid)

def kill_proc():
  chaco_proc.stdin.close()
  chaco_proc.terminate()
  chaco_proc.wait(timeout=0.2)

def write(message):
  chaco_proc.stdin.write(message + "\n")
  chaco_proc.stdin.flush()

t = Timer(args.timelimit, kill_proc)
t.start()
write(str(args.graph))
write("1")                     # Multilevel Kernighan-Lin
write("300")                   # Contraction Limit
write(str(math.log(args.k,2))) # Number of recursive bipartitioning levels
write("1")                     # Recursive Bipartitioning
write("No")                    # Terminate
out, err = chaco_proc.communicate()
t.cancel()
end = time.time()

os.remove("User_Params")

total_time = 2147483647
cut = 2147483647
km1 = 2147483647
soed = 2147483647
imbalance = 1.0
timeout = "no"
failed = "no"

check = False
if chaco_proc.returncode == 0:
  # Extract metrics out of hMetis output
  for line in out.split('\n'):
    s = str(line).strip()
    if ( "nsets = " + str(args.k) in s ):
      check = True
    if (check and "Edge Cuts:" in s):
      cut = int(s.split()[2])
      km1 = cut
    if (check and "partitioning" in s):
      total_time = float(s.split()[1])
    if (check and "Set Size:" in s):
      max_part_size = float(s.split()[3])
      imbalance = max_part_size / math.ceil(float(numNodes)/args.k) - 1.0
elif chaco_proc.returncode == -signal.SIGTERM:
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

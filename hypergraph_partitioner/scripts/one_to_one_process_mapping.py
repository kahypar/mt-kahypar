#!/usr/bin/python3
import subprocess
import argparse
import time
import os
import os.path
from threading import Timer
import signal

###################################
# SETUP ENV
###################################
process_mapping = os.environ.get("PROCESS_MAPPING")
assert (process_mapping != None), "check env.sh"
###################################

parser = argparse.ArgumentParser()
parser.add_argument("graph", type=str)
parser.add_argument("result_file", type=str)

args = parser.parse_args()

result_line = []
with open(args.result_file) as f:
  result_line = f.readline().split(",")
result = { "algorithm": result_line[0],
           "graph": result_line[1],
           "timeout": result_line[2],
           "seed": result_line[3],
           "k": result_line[4],
           "epsilon": result_line[5],
           "num_threads": result_line[6],
           "imbalance": result_line[7],
           "totalPartitionTime": result_line[8],
           "objective": result_line[9],
           "km1": result_line[10],
           "cut": result_line[11],
           "failed": result_line[12].split("\n")[0],
           "totalMappingTime": 2147483647,
           "totalTime": 2147483647,
           "process_mapping": 2147483647,
           "approximation_factor": 2147483647 }

if result["timeout"] == "no" and result["failed"] == "no":
  partition_file = ( os.path.dirname(args.result_file) + "/" + os.path.basename(args.graph) +
                     ".part" + result["k"] + ".epsilon" + result["epsilon"] +
                     ".seed" + result["seed"] + ".partition" )
  process_graph_file = args.graph + ".k" + result["k"]

  if os.path.exists(partition_file) and os.path.exists(process_graph_file):
    # Run One-To-One Process Mapping
    process_mapping_command = [process_mapping,
                              "-h" + args.graph,
                              "-b" + partition_file,
                              "-p" + process_graph_file,
                              "-k" + result["k"],
                              "-s" + result["seed"],
                              "-v0"]
    process_mapping_proc = subprocess.Popen(process_mapping_command,
                            stdout=subprocess.PIPE, universal_newlines=True, preexec_fn=os.setsid)

    def kill_proc():
      os.killpg(os.getpgid(process_mapping_proc.pid), signal.SIGTERM)

    t = Timer(28800, kill_proc)
    t.start()
    out, err = process_mapping_proc.communicate()
    t.cancel()
    end = time.time()

    if process_mapping_proc.returncode == 0:
      for line in out.split('\n'):
        s = str(line).strip()
        if "RESULT" in s:
          result["process_mapping"] = int(s.split(" process_mapping=")[1].split(" ")[0])
          result["approximation_factor"] = float(s.split(" approximation_factor=")[1].split(" ")[0])
          result["totalMappingTime"] = float(s.split(" totalPartitionTime=")[1].split(" ")[0])
          result["totalTime"] = float(result["totalPartitionTime"]) + result["totalMappingTime"]

# CSV format: algorithm,graph,timeout,seed,k,epsilon,num_threads,
# imbalance,totalPartitionTime,totalMappingTime,totalTime,
# objective,km1,cut,process_mapping,approximation_factor,failed
print(result["algorithm"],
      result["graph"],
      result["timeout"],
      result["seed"],
      result["k"],
      result["epsilon"],
      result["num_threads"],
      result["imbalance"],
      result["totalPartitionTime"],
      result["totalMappingTime"],
      result["totalTime"],
      "process_mapping",
      result["km1"],
      result["cut"],
      result["process_mapping"],
      result["approximation_factor"],
      result["failed"],
      sep=",")
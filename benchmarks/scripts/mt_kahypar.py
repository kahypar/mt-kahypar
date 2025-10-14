#!/usr/bin/env python3
import subprocess
import argparse
import sys
import os
import os.path
from threading import Timer
import signal
import shlex
import ntpath
import shutil
import re


_result_values = {
  "timeout": "no",
  "failed": "no",
}
_results_initialized = False


def get_args():
  parser = argparse.ArgumentParser()
  parser.add_argument("graph", type=str)
  parser.add_argument("threads", type=int)
  parser.add_argument("k", type=int)
  parser.add_argument("epsilon", type=float)
  parser.add_argument("seed", type=int)
  parser.add_argument("objective", type=str)
  parser.add_argument("timelimit", type=int)
  parser.add_argument("--partition_folder", type=str, default = "")
  parser.add_argument("--name", type=str, default = "")
  parser.add_argument("--args", type=str, default = "")
  parser.add_argument("--header", type=str, default = "")
  parser.add_argument("--tag", action="store_true")
  return parser.parse_args()

def run_mtkahypar(mt_kahypar, args, default_args, print_fail_msg=True, detect_instance_type=False):
  # Remove --evo marker if present (it's not a Mt-KaHyPar argument)
  cleaned_args = args.args.replace('--evo', '').strip()
  args_list = shlex.split(cleaned_args) if cleaned_args else []

  for arg_key in default_args:
    assert ("--" in arg_key) and not ("=" in arg_key), f"Invalid default argument: {arg_key}"
    if not any(a for a in args_list if (arg_key in a)):
      args_list.append(f"{arg_key}={default_args[arg_key]}")

  if detect_instance_type:
    for arg in args_list:
      assert not ("instance-type" in arg) and not ("input-file-format" in arg)
    if args.graph.endswith(".metis") or args.graph.endswith(".graph"):
      args_list.append("--instance-type=graph")
      args_list.append("--input-file-format=metis")

  #DEBUG PRINT
  #print("DEBUG: Additional Mt-KaHyPar arguments: " + str(args_list), file=sys.stderr)

  # Run Mt-KaHyPar
  cmd = [mt_kahypar,
         "-h" + args.graph,
         "-k" + str(args.k),
         "-e" + str(args.epsilon),
         "--seed=" + str(args.seed),
         "-o" + str(args.objective),
         "-mdirect",
         "--s-num-threads=" + str(args.threads),
         "--verbose=false",
         "--sp-process=true",
         "--show-detailed-timing=true",
         *args_list]
  if args.partition_folder != "":
    cmd.extend(["--write-partition-file=true"])
    cmd.extend(["--partition-output-folder=" + args.partition_folder])
  mt_kahypar_proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True, preexec_fn=os.setsid)

  # handle early interrupt cases where the Mt-KaHyPar process should be killed
  def kill_proc(*args):
    os.killpg(os.getpgid(mt_kahypar_proc.pid), signal.SIGTERM)

  signal.signal(signal.SIGINT, kill_proc)
  signal.signal(signal.SIGTERM, kill_proc)

  t = Timer(args.timelimit, kill_proc)
  t.start()
  out, err = mt_kahypar_proc.communicate()
  t.cancel()

  def _extract_objectives(text: str) -> dict:
      patterns = {
        'km1':       r'^\s*km1\s*=\s*([0-9]+)',
        'cut':       r'^\s*cut\s*=\s*([0-9]+)',
        'soed':      r'^\s*soed\s*=\s*([0-9]+)',
        'imbalance': r'^\s*Imbalance\s*=\s*([0-9.]+)',
        'time':      r'^\s*Partitioning Time\s*=\s*([0-9.]+)\s*s'
      }
      metrics = {}
      
      # DEBUG TEXT
      #print("DEBUG: Mt-KaHyPar output:\n" + text, file=sys.stderr)
      
      for line in text.splitlines():
        for key, pat in patterns.items():
          if key in metrics:
            continue
          m = re.search(pat, line)
          if m:
            if key in ('km1', 'cut', 'soed'):
              metrics[key] = int(m.group(1))
            else:
              metrics[key] = float(m.group(1))
      return metrics
  
  if mt_kahypar_proc.returncode == 0:
    metrics = _extract_objectives(out)
    required = {'km1', 'cut', 'soed', 'imbalance', 'time'}
    
    assert required.issubset(metrics.keys()), "No complete Objectives block found!"
    return metrics, True
  elif mt_kahypar_proc.returncode == -signal.SIGTERM:
    _result_values["timeout"] = "yes"
    return {}, False
  else:
    _result_values["failed"] = "yes"
    if err and print_fail_msg:
      print(err, file=sys.stderr)
  return {}, False
  

def run_mtkahypar_evo(mt_kahypar, args, default_args, print_fail_msg=True, detect_instance_type=False):
  # Remove --evo marker and history-info if present (not Mt-KaHyPar arguments)
  cleaned_args = args.args.replace('--evo', '').strip()
  include_history = False
  if "history-info" in cleaned_args:
    include_history = True
    cleaned_args = cleaned_args.replace('history-info', '').strip()
  args_list = shlex.split(cleaned_args) if cleaned_args else []
  
  #DEBUG PRINT
  #print("DEBUG: Additional Mt-KaHyPar arguments: " + str(args_list), file=sys.stderr)

  for arg_key in default_args:
    assert ("--" in arg_key) and not ("=" in arg_key), f"Invalid default argument: {arg_key}"
    if not any(a for a in args_list if (arg_key in a)):
      args_list.append(f"{arg_key}={default_args[arg_key]}")

  if detect_instance_type:
    for arg in args_list:
      assert not ("instance-type" in arg) and not ("input-file-format" in arg)
    if args.graph.endswith(".metis") or args.graph.endswith(".graph"):
      args_list.append("--instance-type=graph")
      args_list.append("--input-file-format=metis")

  evo_result_file = "/" + ntpath.basename(args.graph) + ".k" + str(args.k) + ".epsilon" + str(args.epsilon) + ".seed" + str(args.seed) + ".timelimit" + str(args.timelimit) + ".csv"
  evo_diff_file = "/" + ntpath.basename(args.graph) + ".k" + str(args.k) + ".epsilon" + str(args.epsilon) + ".seed" + str(args.seed) +  ".timelimit" + str(args.timelimit) + "_diff.csv"

  # Run Mt-KaHyPar
  cmd = [mt_kahypar,
         "-h" + args.graph,
         "-k" + str(args.k),
         "-e" + str(args.epsilon),
         "--seed=" + str(args.seed),
         "-o" + str(args.objective),
         "--s-num-threads=" + str(args.threads),
         "--verbose=false",
         "--sp-process=true",
         "--partition-evolutionary=true",
         "--time-limit=" + str(args.timelimit),
         #"--evo-history-file=" + os.environ.get("EVO_RESULT_FOLDER") + evo_result_file,
         #"--evo-diff-matrix-file=" + os.environ.get("EVO_DIFF_FOLDER") + evo_diff_file,
         *args_list]

  #print("DEBUG: COMMAND: " + " ".join(cmd), file=sys.stderr)

  if args.partition_folder != "":
    cmd.extend(["--write-partition-file=true"])
    cmd.extend(["--partition-output-folder=" + args.partition_folder])
    
  if include_history and os.environ.get("EVO_RESULT_FOLDER") is not None:
    cmd.extend(["--evo-history-file=" + os.environ.get("EVO_RESULT_FOLDER") + evo_result_file])
  if include_history and os.environ.get("EVO_DIFF_FOLDER") is not None:
    cmd.extend(["--evo-diff-matrix-file=" + os.environ.get("EVO_DIFF_FOLDER") + evo_diff_file])
    
  mt_kahypar_proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True, preexec_fn=os.setsid)

  # handle early interrupt cases where the Mt-KaHyPar process should be killed
  def kill_proc(*args):
    os.killpg(os.getpgid(mt_kahypar_proc.pid), signal.SIGTERM)

  signal.signal(signal.SIGINT, kill_proc)
  signal.signal(signal.SIGTERM, kill_proc)

  evo_time_buffer = 10000

  t = Timer(args.timelimit + evo_time_buffer, kill_proc)
  t.start()     
  out, err = mt_kahypar_proc.communicate()
  t.cancel()

  def _extract_objectives(text: str) -> dict:
      patterns = {
        'km1':       r'^\s*km1\s*=\s*([0-9]+)',
        'cut':       r'^\s*cut\s*=\s*([0-9]+)',
        'soed':      r'^\s*soed\s*=\s*([0-9]+)',
        'imbalance': r'^\s*Imbalance\s*=\s*([0-9.]+)',
        'time':      r'^\s*Partitioning Time\s*=\s*([0-9.]+)\s*s'
      }
      metrics = {}
      for line in text.splitlines():
        for key, pat in patterns.items():
          if key in metrics:
            continue
          m = re.search(pat, line)
          if m:
            if key in ('km1', 'cut', 'soed'):
              metrics[key] = int(m.group(1))
            else:
              metrics[key] = float(m.group(1))
      return metrics

  if mt_kahypar_proc.returncode == 0:
    metrics = _extract_objectives(out)
    required = {'km1', 'cut', 'soed', 'imbalance', 'time'}
    assert required.issubset(metrics.keys()), "No complete Objectives block found!"
    return metrics, True
  elif mt_kahypar_proc.returncode == -signal.SIGTERM:
    _result_values["timeout"] = "yes"
    return {}, False
  else:
    _result_values["failed"] = "yes"
    if err and print_fail_msg:
      print(err, file=sys.stderr)
    return {}, False

def set_results(**args):
  global _result_values, _results_initialized
  _result_values.update(args)
  _results_initialized = True

def parse(metrics: dict):
  assert _results_initialized, "set_results must be called before parsing"
  _result_values.update(metrics)
  
def print_result(algorithm, args):
  assert _results_initialized, "set_result_vals must be called before print_result"
  timeout = _result_values["timeout"]
  failed = _result_values["failed"]
  imbalance = _result_values["imbalance"]
  total_time = _result_values["time"]
  km1 = _result_values["km1"]
  cut = _result_values["cut"]
  l = [_result_values]
  del _result_values["timeout"]
  del _result_values["failed"]
  del _result_values["imbalance"]
  del _result_values["time"]
  del _result_values["km1"]
  del _result_values["cut"]
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

  if args.header != "":
    header = ["algorithm", "graph", "timeout", "seed", "k", "epsilon", "num_threads", "imbalance", "totalPartitionTime", "objective", "km1", "cut", "failed"]
    if args.tag:
      header.insert(0, "tag")
    header.extend(_result_values.keys())
    with open(args.header, "w") as header_file:
      header_file.write(",".join(header) + "\n")

  # if args.partition_folder != "":
  #   src_partition_file = args.partition_folder + "/" + ntpath.basename(args.graph) + ".part" + str(args.k) + ".epsilon" + str(args.epsilon) + ".seed" + str(args.seed) + ".KaHyPar"
  #   dst_partition_file = args.partition_folder + "/" + ntpath.basename(args.graph) + ".part" + str(args.k) + ".epsilon" + str(args.epsilon) + ".seed" + str(args.seed) + ".partition"
  #   if os.path.exists(src_partition_file):
  #     shutil.move(src_partition_file, dst_partition_file)


 ###### Main script ######
if __name__ == "__main__":
  algorithm = "Mt-KaHyPar"
  EXECUTABLE = os.environ.get("MT_KAHYPAR")
  assert EXECUTABLE is not None, "Please set the MT_KAHYPAR environment variable to point to the Mt-KaHyPar executable."

  args = get_args()
  if args.name != "":
    algorithm = args.name

  evo = False

  #print(f"DEBUG; CALL WITH ARGS: {args}", file=sys.stderr)

  # determine args evo flag
  if args.args is not None and "evo" in args.args:
    evo = True

  # run Mt-KaHyPar
  result, success = {}, False

  if evo:
    #print("DEBUG: Running Mt-KaHyPar in evolutionary mode.", file=sys.stderr)
    result, success = run_mtkahypar_evo(EXECUTABLE, args, default_args={"--preset-type": "default"}, detect_instance_type=True)
  else:
    result, success = run_mtkahypar(EXECUTABLE, args, default_args={"--preset-type": "default"}, detect_instance_type=True)
  
  set_results()
  if success:
    parse(result)
  print_result(algorithm, args)

  # match(evo):
  #   case True:
  #     result, success = run_mtkahypar_evo(EXECUTABLE, args, default_args={"--preset": "default"}, detect_instance_type=True)
  #   case False:
  #     result, success = run_mtkahypar(EXECUTABLE, args, default_args={"--preset": "default"}, detect_instance_type=True)
  #   case _:
  #     print(f"Unknown mode {args.mode}.", file=sys.stderr)
  #     exit(1)
  # if success:
  #   parse(result)
  # print_result(algorithm, args)
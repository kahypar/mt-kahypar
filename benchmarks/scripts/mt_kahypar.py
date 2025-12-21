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
import threading
import time


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
  parser.add_argument("--history", type=str, default = "")
  parser.add_argument("--iteration_log_file", type=str, default = "")
  parser.add_argument("--args", type=str, default = "")
  parser.add_argument("--header", type=str, default = "")
  parser.add_argument("--tag", action="store_true")
  return parser.parse_args()


def _extract_objectives(text: str) -> dict:
    metrics: dict[str, float | int] = {}

    # First pass: inline key=value anywhere in the text (case-insensitive)
    inline_patterns = {
        'km1':       r'km1\s*=\s*([0-9]+)',
        'cut':       r'cut\s*=\s*([0-9]+)',
        'soed':      r'soed\s*=\s*([0-9]+)',
        'imbalance': r'imbalance\s*=\s*([0-9.]+)',
        'time':      r'(?:totalPartitionTime|totalPartitioningTime)\s*=\s*([0-9.]+)'
    }

    for key, pat in inline_patterns.items():
        if key in metrics:
            continue
        m = re.search(pat, text, flags=re.IGNORECASE)
        if m:
            if key in ('km1', 'cut', 'soed'):
                metrics[key] = int(m.group(1))
            else:
                metrics[key] = float(m.group(1))

    line_patterns = {
        'km1':       r'^\s*KM1\s*=\s*([0-9]+)',
        'cut':       r'^\s*Cut\s*=\s*([0-9]+)',
        'soed':      r'^\s*SOED\s*=\s*([0-9]+)',
        'imbalance': r'^\s*Imbalance\s*=\s*([0-9.]+)%?',
        'time':      r'^\s*(?:Partitioning Time|Total Partitioning Time)\s*=\s*([0-9.]+)\s*s'
    }

    for line in text.splitlines():
        for key, pat in line_patterns.items():
            if key in metrics:
                continue
            m = re.search(pat, line, flags=re.IGNORECASE)
            if m:
                if key in ('km1', 'cut', 'soed'):
                    metrics[key] = int(m.group(1))
                else:
                    metrics[key] = float(m.group(1))

    return metrics


def run_mtkahypar(mt_kahypar, args, default_args, print_fail_msg=True, detect_instance_type=False):
  cleaned_args = args.args.strip() if args.args else ""

  include_history = False
  if cleaned_args and "history-info" in cleaned_args:
    include_history = True
    cleaned_args = cleaned_args.replace("history-info", "").strip()

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

  # evo-style history file for NON-evo runs
  history_path = None
  if include_history and args.history:
    experiment_dir = args.history
    evo_history_dir = os.path.join(experiment_dir, "evo_history")
    os.makedirs(evo_history_dir, exist_ok=True)

    base_filename = ntpath.basename(args.graph)
    unique_suffix = f".k{args.k}.epsilon{args.epsilon}.seed{args.seed}.timelimit{args.timelimit}"
    thread_id = threading.get_ident()
    history_path = os.path.join(
      evo_history_dir,
      f"{base_filename}{unique_suffix}.thread{thread_id}.history.csv"
    )

    with open(history_path, "w", encoding="utf-8") as f:
      f.write(f"Starttime: {int(time.time() * 1000)}\n")

  def _build_cmd(seed: int):
    return [
      mt_kahypar,
      "-h" + args.graph,
      "-k" + str(args.k),
      "-e" + str(args.epsilon),
      "--seed=" + str(seed),
      "-o" + str(args.objective),
      "-mdirect",
      "--s-num-threads=" + str(args.threads),
      "--verbose=false",
      "--sp-process=true",
      "--show-detailed-timing=true",
      *args_list
    ]

  total_budget = float(args.timelimit)
  start = time.monotonic()

  best_metrics = None
  best_km1 = None
  run_idx = 0

  while True:
    elapsed = time.monotonic() - start
    remaining = total_budget - elapsed
    if remaining <= 0:
      break

    seed = args.seed + run_idx
    cmd = _build_cmd(seed)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            universal_newlines=True, preexec_fn=os.setsid)

    def kill_proc(*_):
      try:
        os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
      except Exception:
        pass

    # per-run limit = remaining global budget
    t = Timer(remaining, kill_proc)
    t.start()
    out, err = proc.communicate()
    t.cancel()

    if proc.returncode == 0:
      metrics = _extract_objectives(out)
      required = {'km1', 'cut', 'soed', 'imbalance', 'time'}
      assert required.issubset(metrics.keys()), "No complete Objectives block found!"

      km1 = int(metrics["km1"])
      if best_km1 is None or km1 < best_km1:
        best_km1 = km1
        best_metrics = metrics
        if history_path is not None:
          ts_ms = int(time.time() * 1000)
          with open(history_path, "a", encoding="utf-8") as f:
            # same 3-column format as evo: timestamp, op_type, km1
            f.write(f"{ts_ms}, NORMAL, {km1}\n")
    elif proc.returncode == -signal.SIGTERM:
      # slice timed out; just move on, global loop will stop when budget is gone
      _result_values["timeout"] = "yes"
    else:
      _result_values["failed"] = "yes"
      if err and print_fail_msg:
        print(err, file=sys.stderr)
      # optional: break here instead of continuing if you prefer
    run_idx += 1

  if best_metrics is None:
    return {}, False

  return best_metrics, True


def run_mtkahypar_evo(mt_kahypar, args, default_args, print_fail_msg=True, detect_instance_type=False):
  # Remove --evo marker and history-info if present (not Mt-KaHyPar arguments)
  cleaned_args = args.args.replace('evo_flag', '').strip()
  include_history = False
  iteration_logging = False
  if "history-info" in cleaned_args:
    include_history = True
    cleaned_args = cleaned_args.replace('history-info', '').strip()
  if "iteration-logging" in cleaned_args:
    iteration_logging = True
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

  base_filename = ntpath.basename(args.graph)
  unique_suffix = f".k{args.k}.epsilon{args.epsilon}.seed{args.seed}.timelimit{args.timelimit}"
  thread_id = threading.get_ident()

  evo_result_file = f"/{base_filename}{unique_suffix}.thread{thread_id}.history.csv"
  evo_diff_file = f"/{base_filename}{unique_suffix}.thread{thread_id}.diff.csv"
  evo_iteration_log_file = f"/{base_filename}{unique_suffix}.thread{thread_id}.iteration_log.csv"


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

  #print("DEBUG: args_list: " + str(args_list), file=sys.stderr)

  if args.partition_folder != "":
    cmd.extend(["--write-partition-file=true"])
    cmd.extend(["--partition-output-folder=" + args.partition_folder])

  if include_history:
    # get experiment directory from args
    experiment_dir = args.history
    evo_result_folder = experiment_dir + "/evo_history"
    evo_diff_folder = experiment_dir + "/evo_diff"
    evo_iteration_log_folder = experiment_dir + "/evo_iteration_log"

    # Create directories if they don't exist
    os.makedirs(evo_result_folder, exist_ok=True)
    os.makedirs(evo_diff_folder, exist_ok=True)
    os.makedirs(evo_iteration_log_folder, exist_ok=True)
    
    if evo_result_folder is not None:
      cmd.extend(["--evo-history-file=" + evo_result_folder + evo_result_file])
    if evo_diff_folder is not None:
      cmd.extend(["--evo-diff-matrix-file=" + evo_diff_folder + evo_diff_file])
    if evo_iteration_log_folder is not None and iteration_logging:
      cmd.extend(["--evo-iteration-log-file=" + evo_iteration_log_folder + evo_iteration_log_file])
      
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

  if mt_kahypar_proc.returncode == 0:
    metrics = _extract_objectives(out)
    required = {'km1', 'cut', 'soed', 'imbalance', 'time'}
    assert required.issubset(metrics.keys()), "No complete Objectives block found!"
    return metrics, True
  elif mt_kahypar_proc.returncode == -signal.SIGTERM:
    # DEBUG PRINT
    print(f"DEBUG: Mt-KaHyPar timed out after {args.timelimit} seconds", file=sys.stderr)
    _result_values["timeout"] = "yes"
    return {}, False
  else:
    # DEBUG PRINT
    print(f"DEBUG: Mt-KaHyPar failed with return code {mt_kahypar_proc.returncode}", file=sys.stderr)
    print(f"DEBUG: Mt-KaHyPar stderr output:\n{err}", file=sys.stderr)
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

  user_has_preset = ("--preset-type" in args.args) if args.args else False
  default_args = {} if user_has_preset else {"--preset-type": "quality"}
  
  if evo:
    #print("DEBUG: Running Mt-KaHyPar in evolutionary mode.", file=sys.stderr)
    result, success = run_mtkahypar_evo(EXECUTABLE, args, default_args=default_args, detect_instance_type=True)
  else:
    result, success = run_mtkahypar(EXECUTABLE, args, default_args=default_args, detect_instance_type=True)
  
  set_results()
  if success:
    parse(result)
  else:
    print(f"DEBUG: Mt-KaHyPar failed with result: {result}", file=sys.stderr)
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
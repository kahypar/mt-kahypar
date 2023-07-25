#!/usr/bin/python3
import json
import os
import os.path
import subprocess
import sys
import multiprocessing

mt_kahypar_dir = os.environ.get("PWD") + "/"
executable = mt_kahypar_dir + "build/mt-kahypar/application/MtKaHyPar"
verify_partition_exec = mt_kahypar_dir + "build/tools/VerifyPartition"
generate_grid_graph = mt_kahypar_dir + "build/tools/GridGraphGenerator"
config_dir = mt_kahypar_dir + "config/"
integration_test_json_file = mt_kahypar_dir + "tests/end_to_end/integration_tests.json"
num_threads = multiprocessing.cpu_count()


partitioners = { "Mt-KaHyPar-D":       { "config":  config_dir + "default_preset.ini",
                                         "mode": "direct" },
                 "Mt-KaHyPar-Q":       { "config":  config_dir + "quality_preset.ini",
                                         "mode": "direct" },
                 "Mt-KaHyPar-Q-F":     { "config":  config_dir + "highest_quality_preset.ini",
                                         "mode": "direct" },
                 "Mt-KaHyPar-D-Graph": { "config":  config_dir + "default_preset.ini",
                                         "mode": "direct" },
                 "Mt-KaHyPar-Q-Graph": {"config":  config_dir + "quality_preset.ini",
                                         "mode": "direct" },
                 "Mt-KaHyPar-Q-F-Graph": {"config":  config_dir + "highest_quality_preset.ini",
                                         "mode": "direct" },
                 "Mt-KaHyPar-SDet":    { "config":  config_dir + "deterministic_preset.ini",
                                         "mode": "direct" },
                 "Mt-KaHyPar-D-RB":    { "config":  config_dir + "default_preset.ini",
                                         "mode": "rb" },
                 "Mt-KaHyPar-Q-RB":    { "config":  config_dir + "quality_preset.ini",
                                         "mode": "rb" },
                 "Mt-KaHyPar-Q-F-RB":  { "config":  config_dir + "highest_quality_preset.ini",
                                         "mode": "rb" },
                 "Mt-KaHyPar-D-Deep":  { "config":  config_dir + "default_preset.ini",
                                         "mode": "deep" },
                 "Mt-KaHyPar-Q-Deep":  { "config":  config_dir + "quality_preset.ini",
                                         "mode": "deep" },
                 "Mt-KaHyPar-Q-F-Deep": { "config":  config_dir + "highest_quality_preset.ini",
                                          "mode": "deep" } }

def bold(msg):
  return "\033[1m" + msg + "\033[0m"

def print_error(msg):
  print("\033[1;91m[ERROR]\033[0m " + bold(msg))

def print_success(msg):
  print("\033[1;92m[SUCCESS]\033[0m " + bold(msg))

def print_config(instance, k, epsilon):
  print(bold("Instance = " + instance + ", k = " + str(k) + ", Epsilon = " + str(epsilon)))

def command(test, instance, k, epsilon, target_graph):
  partitioner = partitioners[test["partitioner"]]
  config = config_dir + test["config"] if "config" in test else partitioner["config"]
  parameters = test["parameters"] if "parameters" in test else []
  cmd = [ executable,
         "-h" + instance,
         "-p" + config,
         "-k" + str(k),
         "-e" + str(epsilon),
         "-t" + str(num_threads),
         "-m" + partitioner["mode"],
         "--seed=1",
         "--show-detailed-timings=false",
         "--sp-process=true",
         "--write-partition-file=true"] + parameters
  if target_graph == "":
    cmd = cmd + ["-okm1"]
  else:
    cmd = cmd + ["-osteiner_tree", "-g" + target_graph]
  return cmd

def grep_result(out):
  is_mapping = False
  for line in out.split('\n'):
    s = str(line).strip()
    if "RESULT" in s:
      km1 = int(s.split(" km1=")[1].split(" ")[0])
      soed = int(s.split(" soed=")[1].split(" ")[0])
      cut = int(s.split(" cut=")[1].split(" ")[0])
      if " steiner_tree=" in s:
        steiner_tree = int(s.split(" steiner_tree=")[1].split(" ")[0])
        is_mapping = True
      total_time = float(s.split(" totalPartitionTime=")[1].split(" ")[0])
      imbalance = float(s.split(" imbalance=")[1].split(" ")[0])
  res = "Total Time = " + str(total_time) + \
        " Imbalance = " + str(imbalance) + \
        " km1 = " + str(km1) + \
        " cut = " + str(cut) + \
        " soed = " + str(soed)
  if is_mapping:
    res = res + " steiner_tree = " + str(steiner_tree)
  return res

def partition_file_str(instance, k, epsilon):
  return instance + ".part" + str(k) + ".epsilon" + str(epsilon) + ".seed1.KaHyPar"

def verify_partition(instance, k, epsilon):
  partition_file = partition_file_str(instance, k, epsilon)
  proc = subprocess.Popen([verify_partition_exec,
                           "-h" + instance,
                           "-b" + partition_file,
                           "-k" + str(k),
                           "-e" + str(epsilon)],
                           stdout=subprocess.PIPE, universal_newlines=True)
  out, err = proc.communicate()
  os.remove(partition_file)

  if proc.returncode == 0:
    print_success("Partition State: VALID")
  else:
    print_error("Partition State: INVALID")
    print(out)
    sys.exit(-1)

def generate_target_graph(instance, k):
  n = int(k / 2)
  m = 2
  tmp_k = n * m
  proc = subprocess.Popen([generate_grid_graph,
                           "-o" + instance,
                           "--n=" + str(n),
                           "--m=" + str(m),
                           "--max-weight=10"],
                           stdout=subprocess.PIPE, universal_newlines=True)
  out, err = proc.communicate()

  if proc.returncode == 0:
    return instance + ".k" + str(tmp_k)
  else:
    return ""

def run_integration_test(cmd, instance, k, epsilon):
  proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
  out, err = proc.communicate()

  if proc.returncode == 0:
    verify_partition(instance, k, epsilon)
    print_success(grep_result(out))
  else:
    print_error("Partitioner terminates with non-zero exit code (Exit Code = " + str(proc.returncode) + ")")
    print(out)
    sys.exit(-1)


with open(integration_test_json_file) as integration_test_file:
  integration_tests = json.load(integration_test_file)
  for experiment in integration_tests["integration_tests"]:
    print(bold(experiment["name"]))
    print("".rjust(len(experiment["name"]), "-"))
    epsilon = integration_tests["epsilon"]
    for instance in experiment["instances"]:
      absoulute_instance_path = mt_kahypar_dir + instance
      for k in integration_tests["k"]:
        for test in experiment["tests"]:
          target_graph = ""
          if "mapping" in test:
            if k % 2 != 0:
              continue
            target_graph = generate_target_graph(absoulute_instance_path, k)
          print_config(instance, k, epsilon)
          cmd = command(test, absoulute_instance_path, k, epsilon, target_graph)
          print(' '.join(cmd))
          run_integration_test(cmd, absoulute_instance_path, k, epsilon)
          if target_graph != "" and os.path.exists(target_graph):
            os.remove(target_graph)
          print()




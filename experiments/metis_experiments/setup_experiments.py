#!/usr/bin/python3
import json
import argparse
import datetime
import os
import os.path
import ntpath
import shutil
import re

def intersection(lst1, lst2):
  lst3 = [value for value in lst1 if value in lst2]
  return lst3

partitioner_script_folder = os.environ.get("PARTITIONER_SCRIPT_FOLDER")
assert (partitioner_script_folder != None), "check env.sh"

serial_partitioner = [ "hMetis-R", "hMetis-K", "PaToH-S", "PaToH-D", "PaToH-Q",
                       "KaHyPar-CA", "KaHyPar-K", "KaHyPar-R", "Mondriaan", "Hype",
                       "KaFFPa-Fast", "KaFFPa-FastS", "KaFFPa-Eco", "KaFFPa-EcoS",
                       "KaFFPa-Strong", "KaFFPa-StrongS", "Metis-R", "Metis-K",
                       "Chaco-R", "Chaco-K", "Scotch" ]
parallel_partitioner = [ "Parkway", "Zoltan", "MT-KaHyPar-D", "MT-KaHyPar-Q", "MT-KaHyPar-D-F", "MT-KaHyPar-Q-F",
                         "MT-KaHyPar-Graph-D", "MT-KaHyPar-Graph-Q", "MT-KaHIP", "MT-Metis", "KaMinPar",
                         "ParHIP", "ParMetis", "PT-Scotch", "BiPart" ]

partitioner_mapping = { "hMetis-R": "hmetis_rb",
                        "hMetis-K": "hmetis_k",
                        "PaToH-S": "patoh_s",
                        "PaToH-D": "patoh_d",
                        "PaToH-Q": "patoh_q",
                        "KaHyPar-CA": "kahypar_ca",
                        "KaHyPar-K": "kahypar_k",
                        "KaHyPar-R": "kahypar_r",
                        "Mondriaan": "mondriaan",
                        "Hype": "hype",
                        "Parkway": "parkway",
                        "Zoltan": "zoltan",
                        "BiPart": "bipart",
                        "MT-KaHyPar": "mt_kahypar",
                        "MT-KaHyPar-D": "mt_kahypar_d",
                        "MT-KaHyPar-Q": "mt_kahypar_q",
                        "MT-KaHyPar-D-F": "mt_kahypar_d_f",
                        "MT-KaHyPar-Q-F": "mt_kahypar_q_f",
                        "MT-KaHyPar-Graph-D": "mt_kahypar_graph_d",
                        "MT-KaHyPar-Graph-Q": "mt_kahypar_graph_q",
                        "MT-KaHIP": "mt_kahip",
                        "MT-Metis": "mt_metis",
                        "KaMinPar": "kaminpar",
                        "Metis-R": "metis_rb",
                        "Metis-K": "metis_k",
                        "Chaco-R": "chaco_rb",
                        "Chaco-K": "chaco_k",
                        "Scotch": "scotch",
                        "PT-Scotch": "pt_scotch",
                        "KaFFPa-Fast": "kaffpa_fast",
                        "KaFFPa-FastS": "kaffpa_fastsocial",
                        "KaFFPa-Eco": "kaffpa_eco",
                        "KaFFPa-EcoS": "kaffpa_ecosocial",
                        "KaFFPa-Strong": "kaffpa_strong",
                        "KaFFPa-StrongS": "kaffpa_strongsocial",
                        "ParHIP": "parhip",
                        "ParMetis": "parmetis",
                        "multi_dim_mt_kahypar": "mt_kahypar_md"}

format_mapping = { "mc_metis": "mc_metis_instance_folder",
                   "multi_dim_mt_kahypar": "multi_dim_mt_kahypar_instance_folder", 
                   "hMetis-R": "hmetis_instance_folder",
                   "hMetis-K": "hmetis_instance_folder",
                   "PaToH-S": "patoh_instance_folder",
                   "PaToH-D": "patoh_instance_folder",
                   "PaToH-Q": "patoh_instance_folder",
                   "KaHyPar-CA": "hmetis_instance_folder",
                   "KaHyPar-K": "hmetis_instance_folder",
                   "KaHyPar-R": "hmetis_instance_folder",
                   "Mondriaan": "hmetis_instance_folder",
                   "Hype": "hmetis_instance_folder",
                   "Parkway": "hmetis_instance_folder",
                   "Zoltan": "zoltan_instance_folder",
                   "BiPart": "hmetis_instance_folder",
                   "MT-KaHyPar-D": "hmetis_instance_folder",
                   "MT-KaHyPar-Q": "hmetis_instance_folder",
                   "MT-KaHyPar-D-F": "hmetis_instance_folder",
                   "MT-KaHyPar-Q-F": "hmetis_instance_folder",
                   "MT-KaHyPar-Graph-D": "graph_instance_folder",
                   "MT-KaHyPar-Graph-Q": "graph_instance_folder",
                   "MT-KaHIP": "graph_instance_folder",
                   "MT-Metis": "graph_instance_folder",
                   "KaMinPar": "graph_instance_folder",
                   "Metis-R": "graph_instance_folder",
                   "Metis-K": "graph_instance_folder",
                   "Chaco-R": "graph_instance_folder",
                   "Chaco-K": "graph_instance_folder",
                   "Scotch": "scotch_instance_folder",
                   "PT-Scotch": "scotch_instance_folder",
                   "KaFFPa-Fast": "graph_instance_folder",
                   "KaFFPa-FastS": "graph_instance_folder",
                   "KaFFPa-Eco": "graph_instance_folder",
                   "KaFFPa-EcoS": "graph_instance_folder",
                   "KaFFPa-Strong": "graph_instance_folder",
                   "KaFFPa-StrongS": "graph_instance_folder",
                   "ParHIP": "graph_instance_folder",
                   "ParMetis": "graph_instance_folder" }

def get_all_hypergraph_instances(dir):
  return [dir + "/" + hg for hg in os.listdir(dir) if hg.endswith('.hgr')]

def get_all_mondriaan_instances(dir):
  return [dir + "/" + mondriaan_hg for mondriaan_hg in os.listdir(dir) if mondriaan_hg.endswith('.mondriaan.mtx')]

def get_all_zoltan_instances(dir):
  return [dir + "/" + zoltan_hg for zoltan_hg in os.listdir(dir) if zoltan_hg.endswith('.zoltan.hg')]

def get_all_graph_instances(dir):
  return [dir + "/" + graph for graph in os.listdir(dir) if graph.endswith('.graph')]

def get_all_scotch_instances(dir):
  return [dir + "/" + graph for graph in os.listdir(dir) if graph.endswith('.scotch')]

def get_all_benchmark_instances(partitioner, config):
  config_instance_type = format_mapping[partitioner]
  instance_dir = config[config_instance_type]
  if partitioner == "Mondriaan":
    return get_all_mondriaan_instances(instance_dir)
  elif config_instance_type == "hmetis_instance_folder" or config_instance_type == "patoh_instance_folder":
    return get_all_hypergraph_instances(instance_dir)
  elif config_instance_type == "zoltan_instance_folder":
    return get_all_zoltan_instances(instance_dir)
  elif config_instance_type == "graph_instance_folder":
    return get_all_graph_instances(instance_dir)
  elif config_instance_type == "scotch_instance_folder":
    return get_all_scotch_instances(instance_dir)
  elif config_instance_type == "mc_metis_instance_folder":
    return get_all_graph_instances(instance_dir)
  elif config_instance_type == "multi_dim_mt_kahypar_instance_folder":
    return get_all_hypergraph_instances(instance_dir)

def serial_partitioner_call(partitioner, instance, k, epsilon, seed, objective, timelimit, config_file, algorithm_name, args):
  call = partitioner_script_folder + "/" + partitioner_mapping[partitioner] + ".py " + " " + instance + " " + str(k) + " " +"\"" + str(epsilon) +"\"" " " + str(seed) + " " + str(objective) + " " + str(timelimit)
  if config_file != "":
    call = call + " --config " + config_file
  if algorithm_name != "":
    call = call + " --name " + algorithm_name
  if args != None:
    assert "'" not in args
    call = call + f" --args ' {args}'"
  return call

def parallel_partitioner_call(partitioner, instance, threads, k, epsilon, seed, objective, timelimit, config_file, algorithm_name, args):
  call = partitioner_script_folder + "/" + partitioner_mapping[partitioner] + ".py " + instance + " " + str(threads) + " " + str(k) + " " + "\"" + str(epsilon) + "\"" + " " + str(seed) + " " + str(objective) + " " + str(timelimit)
  if config_file != "":
    call = call + " --config " + config_file
  if algorithm_name != "":
    call = call + " --name " + algorithm_name
  if args != None:
    assert "'" not in args
    call = call + f" --args ' {args}'"
  return call

def partitioner_dump(result_dir, instance, threads, k, seed):
  return os.path.abspath(result_dir) + "/" + ntpath.basename(instance) + "." + str(threads) + "." + str(k) + "." + str(seed) + ".results"

parser = argparse.ArgumentParser()
parser.add_argument("experiment", type=str)

args = parser.parse_args()

with open(args.experiment) as json_experiment:
    config = json.load(json_experiment)

    now = datetime.datetime.now()
    experiment_dir = str(now.year) + "-" + str(now.month) + "-" + str(now.day) + "_" + config["name"]
    workload_file = experiment_dir + "/workload.txt"
    shutil.rmtree(experiment_dir, ignore_errors=True)
    os.makedirs(experiment_dir, exist_ok=True)

    epsilon = config["epsilon"]
    objective = config["objective"]
    timelimit = config["timelimit"]
    write_partition_file = config["write_partition_file"] if "write_partition_file" in config else False

    # Setup experiments
    for partitioner_config in config["config"]:
      partitioner = partitioner_config["partitioner"]
      algorithm_file = partitioner
      if "name" in partitioner_config:
        algorithm_file = partitioner_config["name"]
      algorithm_file = '_'.join(list(map(lambda x: x.lower(), re.split(' |-', algorithm_file))))
      result_dir = experiment_dir + "/" + algorithm_file + "_results"
      os.makedirs(result_dir, exist_ok=True)

    for seed in config["seeds"]:
      for partitioner_config in config["config"]:
        partitioner = partitioner_config["partitioner"]
        algorithm_file = partitioner
        if "name" in partitioner_config:
          algorithm_file = partitioner_config["name"]
        algorithm_file = '_'.join(list(map(lambda x: x.lower(), re.split(' |-', algorithm_file))))
        result_dir = experiment_dir + "/" + algorithm_file + "_results"

        is_serial_partitioner = partitioner in serial_partitioner
        config_file = ""
        if "config_file" in partitioner_config:
          config_file = partitioner_config["config_file"]
        algorithm_name = '"' + partitioner + '"'
        if "name" in partitioner_config:
          algorithm_name = '"' + partitioner_config["name"] + '"'
        args = None
        if "args" in partitioner_config:
          args = partitioner_config["args"]

        for instance in get_all_benchmark_instances(partitioner, config):
          for k in config["k"]:
              partitioner_calls = []
              for threads in config["threads"]:
                if is_serial_partitioner and threads > 1 and len(config["threads"]) > 1:
                  continue
                if is_serial_partitioner:
                  partitioner_call = serial_partitioner_call(partitioner, instance, k, epsilon, seed, objective, timelimit, config_file, algorithm_name, args)
                else:
                  partitioner_call = parallel_partitioner_call(partitioner, instance, threads, k, epsilon, seed, objective, timelimit, config_file, algorithm_name, args)
                if write_partition_file:
                  partitioner_call = partitioner_call + " --partition_folder=" + os.path.abspath(result_dir)
                partitioner_call = partitioner_call + " >> " + partitioner_dump(result_dir, instance, threads, k, seed)
                partitioner_calls.extend([partitioner_call])

              # Write partitioner calls to workload file
              with open(experiment_dir + "/" + algorithm_file + "_workload.txt", "w") as partitioner_workload_file:
                partitioner_workload_file.write("\n".join(partitioner_calls))
                partitioner_workload_file.write("\n")

              with open(workload_file, "a") as global_workload_file:
                global_workload_file.write("\n".join(partitioner_calls))
                global_workload_file.write("\n")



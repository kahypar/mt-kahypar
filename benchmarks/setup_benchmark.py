#!/usr/bin/python3
import json
import argparse
import datetime
import os
import os.path
import ntpath
import shutil
import re

from partitioner_mapping import partitioner_mapping

partitioner_script_folder = os.environ.get("PARTITIONER_SCRIPT_FOLDER")

print("Using partitioner scripts from: " + str(partitioner_script_folder))

assert (partitioner_script_folder != None), "check env.sh"

def get_all_hypergraph_instances(dir):
  return [dir + "/" + hg for hg in os.listdir(dir) if hg.endswith('.hgr') or hg.endswith('.hmetis') or hg.endswith('.graph')]

def partitioner_header(result_dir):
  return str(os.path.abspath(result_dir)).removesuffix("_results") + ".header.csv"

def serial_partitioner_call(partitioner, instance, k, epsilon, seed, objective, timelimit):
  return (
    partitioner_script_folder + "/" + partitioner_mapping[partitioner].script + ".py " + instance
    + " " + str(k) + " " + str(epsilon) + " " + str(seed) + " " + str(objective) + " " + str(timelimit)
  )

def parallel_partitioner_call(partitioner, instance, threads, k, epsilon, seed, objective, timelimit):
  return (
    partitioner_script_folder + "/" + partitioner_mapping[partitioner].script + ".py " + instance
    + " " + str(threads) + " " + str(k) + " " + str(epsilon) + " " + str(seed) + " " + str(objective) + " " + str(timelimit)
  )

def get_all_benchmark_instances(partitioner, config):
  input_format_list = partitioner_mapping[partitioner].format
  for format_type in input_format_list:
    input_format = format_type + "_instance_folder"
    if input_format in config:
      assert "instances" not in config
      instance_dir = config[input_format]
      return {instance: None for instance in get_all_hypergraph_instances(instance_dir)}

    elif "instances" in config:
      # more general case where multiple directories and tags can be defined
      result = {}
      dir_list = [(entry["path"], entry["type"], entry.get("tag")) for entry in config["instances"]]
      for (instance_dir, curr_format, instance_tag) in dir_list:
        assert curr_format in ["hmetis", "patoh", "zoltan", "graph", "metis", "scotch"], f"invalid instance type: {curr_format}"
        if curr_format == format_type:
          input_format = format_type + "_instance_folder"
          tmp = {instance: instance_tag for instance in get_all_hypergraph_instances(instance_dir)}
          intersection = {ntpath.basename(p) for p in result} & {ntpath.basename(p) for p in tmp}
          assert len(intersection) == 0, f"instance appears in multiple folders: {intersection}"
          result.update({instance: instance_tag for instance in tmp})
      if len(result) > 0:
        assert all(tag is not None for tag in result.values()) or all(tag is None for tag in result.values()), "Inconsistent instance tags!"
        return result

  assert False, f"No instances found for: {partitioner}"

def partitioner_call(is_serial, partitioner, instance, threads, k, epsilon, seed, objective, timelimit, config_file, algorithm_name, args, header, tag, experiment_dir):
  if is_serial:
    call = serial_partitioner_call(partitioner, instance, k, epsilon, seed, objective, timelimit)
  else:
    call = parallel_partitioner_call(partitioner, instance, threads, k, epsilon, seed, objective, timelimit)
  if config_file != "":
    call += " --config " + config_file
  if algorithm_name != "":
    call += " --name " + algorithm_name
  if args is not None:
    assert "'" not in args
    if "history-info" in args:
       # pass experiment directory for evo history info
       call += f" --history " + experiment_dir
    call += f" --args '{args}'"
  if header is not None:
    call += f" --header '{header}'"
    if tag is not None:
      call += " --tag"
  if tag is not None:
    call += f' | {{ line=$(cat); echo "{tag},$line"; }}'  # bash snippet which prepends to stdin
  return call

def partitioner_dump(result_dir, instance, threads, k, seed, timelimit):
  return os.path.abspath(result_dir) + "/" + ntpath.basename(instance) + "." + str(threads) + "." + str(k) + "." + str(seed) + "." + str(timelimit) + ".results"

### MAIN SCRIPT ###
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("benchmark", type=str)
    parser.add_argument("-f", "--force", action="store_true")

    args = parser.parse_args()

    with open(args.benchmark, 'r') as f:
        config = json.load(f)
    
    now = datetime.datetime.now()
    experiment_dir = "benchmarks/" + str(now.year) + "-" + str(now.month) + "-" + str(now.day) + "_" + config["name"]
    workload_file = experiment_dir + "/workload.txt"
    if args.force:
        shutil.rmtree(experiment_dir, ignore_errors=True)
        os.makedirs(experiment_dir, exist_ok=True)
    else:
        try:
            os.makedirs(experiment_dir, exist_ok=False)
        except OSError:
            print("Experiment directory already exists! Call with -f to delete old directory")
            exit(1)

    epsilon = config["epsilon"]
    objective = config["objective"]
    timelimit = config["timelimit"]
    write_partition_file = config["write_partition_file"] if "write_partition_file" in config else False
    dynamic_header = config["dynamic_header"] if "dynamic_header" in config else True
    
    repetition_compare_benchmark = False
    if config["repetitions"] is not None and len(config["repetitions"]) > 1:
        repetition_compare_benchmark = True
        assert len(config["repetitions"]) == len(config["timelimit"]), "If repetitions are specified, the number of repetitions must be equal to the number of threads"

    # Setup experiments
    try:
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

            is_serial_partitioner = not partitioner_mapping[partitioner].parallel
            config_file = ""
            if "config_file" in partitioner_config:
                config_file = partitioner_config["config_file"]
            algorithm_name = '"' + partitioner + '"'
            if "name" in partitioner_config:
                algorithm_name = '"' + partitioner_config["name"] + '"'
            args = None
            if "args" in partitioner_config:
                args = partitioner_config["args"]
            header = None
            if dynamic_header and partitioner_mapping[partitioner].dynamic_header:
                header = partitioner_header(result_dir)

            partitioner_calls = []
            for instance, tag in get_all_benchmark_instances(partitioner, config).items():
                for k in config["k"]:
                    for threads in config["threads"]:  
                      if not repetition_compare_benchmark:
                        if is_serial_partitioner and threads > 1 and len(config["threads"]) > 1:
                              continue
                        timelimit = config["timelimit"][0]    
                        call = partitioner_call(is_serial_partitioner, partitioner, instance, threads, k, epsilon, seed, objective, timelimit, config_file, algorithm_name, args, header, tag, experiment_dir)
                        header = None
                        if write_partition_file:
                            call += " --partition_folder=" + os.path.abspath(result_dir)
                        call += " >> " + partitioner_dump(result_dir, instance, threads, k, seed, timelimit)
                        partitioner_calls.append(call)
                      else:
                        # execute timelimit[i] with repetitions[i]
                        for i, repetitions in enumerate(config["repetitions"]):
                            timelimit = config["timelimit"][i]
                            if is_serial_partitioner and threads > 1 and len(config["threads"]) > 1:
                                continue
                            for r in range(repetitions):
                                real_seed = seed + r
                                call = partitioner_call(is_serial_partitioner, partitioner, instance, threads, k, epsilon, real_seed, objective, timelimit, config_file, algorithm_name, args, header, tag, experiment_dir)
                                header = None
                                if write_partition_file:
                                    call += " --partition_folder=" + os.path.abspath(result_dir)
                                call += " >> " + partitioner_dump(result_dir, instance, threads, k, real_seed, timelimit)
                                partitioner_calls.append(call)

            # Write partitioner calls to workload file
            with open(experiment_dir + "/" + algorithm_file + "_workload.txt", "a") as partitioner_workload_file:
                partitioner_workload_file.write("\n".join(partitioner_calls))
                partitioner_workload_file.write("\n")

            with open(workload_file, "a") as global_workload_file:
                global_workload_file.write("\n".join(partitioner_calls))
                global_workload_file.write("\n")

    except AssertionError as e:
        shutil.rmtree(experiment_dir, ignore_errors=True)
        raise e
    except FileNotFoundError as e:
        shutil.rmtree(experiment_dir, ignore_errors=True)
        raise e
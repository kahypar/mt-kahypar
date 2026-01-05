#!/usr/bin/python3
import json
import argparse
import datetime
import os
import os.path
import re


# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = "\033[1;92m" + fill * filledLength + "\033[0m" + ' ' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix} ({iteration}/{total})', end = printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()

parser = argparse.ArgumentParser()
parser.add_argument("benchmark", type=str)

args = parser.parse_args()

with open(args.benchmark) as json_experiment:
    config = json.load(json_experiment)

    now = datetime.datetime.now()
    experiment_dir = "benchmarks/" +  str(now.year) + "-" + str(now.month) + "-" + str(now.day) + "_" + config["name"]
    workload_file = experiment_dir + "/workload.txt"
    num_lines = sum(1 for line in open(workload_file))

    with open(workload_file) as workload:
      lines = workload.readlines()
    lines = [line.rstrip() for line in lines]

    for partitioner_config in config['config']:
      partitioner = partitioner_config["partitioner"]
      algorithm_name = partitioner
      if "name" in partitioner_config:
        algorithm_name = partitioner_config["name"]
      algorithm_name = "_".join(map(lambda x: x.lower(), re.split(r"[ \-]", algorithm_name)))
      os.system("rm " + experiment_dir + "/" + algorithm_name + "_results/*")

    # GNU parallel calls
    # get phyical CPU core count
    cpu_cores = os.cpu_count() // 2
    max_threads_per_call = max(config["threads"])

    # allow overriding the number of parallel jobs via config["jobs"]
    jobs_override = config.get("jobs", None)
    if jobs_override is not None:
        final_cores = int(jobs_override)
    else:
        final_cores = min(cpu_cores // max_threads_per_call, 42)  # limit to 42 parallel jobs to avoid overloading the system

    print(f"Executing benchmark with {final_cores} physical CPU cores and max {max_threads_per_call} threads per partitioner call.")
    os.system(f"parallel --jobs {final_cores} --memfree 5G  :::: {workload_file}")
    

    # previous sequential logic
    # i = 0
    # for partitioner_call in lines:
    #   #print("BEFORE PARTITIONER CALL")
    #   printProgressBar(i, num_lines, prefix = "Progress:", suffix = "Completed")
    #   os.system(partitioner_call)
    #   #print("AFTER PARTITIONER CALL")
    #   i = i + 1
    # printProgressBar(num_lines, num_lines, prefix = "Progress:", suffix = "Completed")

    for partitioner_config in config['config']:
      partitioner = partitioner_config["partitioner"]
      algorithm_name = partitioner
      if "name" in partitioner_config:
        algorithm_name = partitioner_config["name"]
      algorithm_name = '_'.join(list(map(lambda x: x.lower(), re.split(' |-', algorithm_name))))
      result_file = experiment_dir + "/" + algorithm_name + ".csv"
      header_file = experiment_dir + "/" + algorithm_name + ".header.csv"
      if os.path.exists(result_file):
        os.remove(result_file)
      if os.path.exists(header_file):
        os.system("cat \"" + header_file + "\" >> " + result_file)
      else:
        os.system("echo 'algorithm,graph,timeout,seed,k,epsilon,num_threads,imbalance,totalPartitionTime,objective,km1,cut,failed, soed' >> " + result_file)
      os.system("cat " + experiment_dir + "/" + algorithm_name + "_results/* >> " + result_file)



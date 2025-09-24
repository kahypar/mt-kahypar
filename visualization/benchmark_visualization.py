from datetime import datetime
from os import times
from matplotlib.ticker import FuncFormatter, ScalarFormatter, FormatStrFormatter
import matplotlib.pyplot as plt
import pandas as pd
import math
import numpy as np
import os
from enum import Enum

BASELINE_TIME = .3
normalize_time = lambda n: n / BASELINE_TIME

class Mode(Enum):
    DEFAULT = "default"
    EVO = "evo"

def parse_benchamark_input(input_file, mode: Mode):
    lowest_value = math.inf
    seed_index = 0
    with open(input_file, 'r', encoding='ascii') as file:
        lines = [line.strip() for line in file.readlines()]
        seed_count = lines.count('###') + 1
        all_seeds_results = [[] for _ in range(seed_count)]
        start_time = datetime.strptime(lines[1][:-3], "%H:%M:%S.%f")
        last_time = datetime.strptime(lines[1][:-3], "%H:%M:%S.%f")
        print(f"Start Time: {start_time} run {seed_index}")
        seed_results = [[], []]
        max_results = 0
        results = 0
        all_seeds_results = [[] for _ in range(seed_count)]
        lowest_value = math.inf
        for line in lines[2:]:
            if(line == '###'):
                lowest_value = math.inf
                seed_index += 1
                all_seeds_results[seed_index - 1] = seed_results
                seed_results = [[], []]
                results = 0
                continue
            elif ',' not in line:
                start_time = datetime.strptime(line[:-3], "%H:%M:%S.%f")
                last_time = datetime.strptime(line[:-3], "%H:%M:%S.%f")
                print('-'*20)
                print(f"Start Time: {start_time} run {seed_index}")
                continue
            timestamp, method, traget_value = line.split(', ')
            traget_value = int(traget_value)
            time = datetime.strptime(timestamp[:-3], "%H:%M:%S.%f")
            delta = None
            if mode == Mode.DEFAULT:
                delta = (time - start_time).total_seconds()
            elif mode == Mode.EVO:
                delta = (time - last_time).total_seconds()
                last_time = time
            else:
                print("Not implemented")
            if(traget_value < lowest_value):
                lowest_value = traget_value
                seed_results[0] +=  [normalize_time(delta)]
                seed_results[1] += [traget_value]
                print(f"Time: {delta}, Method: {method}, Target Value: {traget_value}")
                results += 1
                max_results = max(results, max_results)
        all_seeds_results[-1] = seed_results
        first_result = [(max(res[0][0] for res in all_seeds_results), np.mean([res[1][0] for res in all_seeds_results]))]
        result_array = first_result + [(math.inf,0) for i in range(sum(len(j[0])for j in all_seeds_results)-seed_count)] 
        print('Done')

        seed_min_indices = [1 if i[0][0] != first_result[0][0] else 1 for i in all_seeds_results]
        for id in range(1, len(result_array)):
            # find the next minimum timestamp
            min_timestamp = math.inf
            min_seed_index = -1
            for seed_i in range(seed_count):
                if seed_min_indices[seed_i] < len(all_seeds_results[seed_i][0]):
                    if all_seeds_results[seed_i][0][seed_min_indices[seed_i]] < min_timestamp:
                        min_timestamp = all_seeds_results[seed_i][0][seed_min_indices[seed_i]]
                        min_seed_index = seed_i

            if min_seed_index != -1:
                seed_min_indices[min_seed_index] += 1
                average = np.mean([all_seeds_results[i][1][seed_min_indices[i]-1] if seed_min_indices[i] > 0 and seed_min_indices[i]-1 < len(all_seeds_results[i][1]) else all_seeds_results[i][1][-1] for i in range(seed_count)])
                result_array[id] = (min_timestamp, average)
        print(result_array)

    return result_array
    
def create_plot(datasets, title, xlabel, ylabel, output_file, scientific=False):
    """
    datasets: list of tuples [(data, label), (data, label), ...]
    where data is [(time, value), (time, value), ...]
    """
    plt.figure(figsize=(10, 6))
    
    if isinstance(datasets[0], tuple) and isinstance(datasets[0][0], (int, float)):
        datasets = [(datasets, "Data")]
    
    for data, label in datasets:
        times = [d[0] for d in data]
        values = [d[1] for d in data]
        plt.plot(times, values, marker='o', label=label)
    
    plt.xscale('log')
    # plt.yscale('log')
    
    ax = plt.gca()
    
    if not scientific:
        # Force plain decimal format for both axes
        for axis in [ax.xaxis, ax.yaxis]:
            formatter = ScalarFormatter()
            formatter.set_scientific(False)
            axis.set_major_formatter(formatter)
            axis.set_minor_formatter(formatter)
    else:
        # Use default scientific notation for log scale
        pass # Matplotlib's default for log scale is scientific

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.savefig(output_file)
    plt.close()

# Multiple datasets:
script_dir = os.path.dirname(os.path.abspath(__file__))

# Use relative paths for input files (assuming they're in the same directory)
result1 = parse_benchamark_input('output.csv', mode=Mode.DEFAULT)
result2 = parse_benchamark_input('output_evo.csv', mode=Mode.EVO)
datasets = [(result1, "Evolutionary"), (result2, "Baseline")]

# Save output file in the same directory as the script
output_path = os.path.join(script_dir, "comparison.png")
create_plot(datasets, "Benchmark Comparison", "Time (t_n)", "Value", output_path, scientific=False)

print(f"Plot saved to: {output_path}")
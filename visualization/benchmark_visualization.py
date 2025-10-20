from datetime import datetime
from os import times
from matplotlib.pylab import seed
from matplotlib.ticker import FuncFormatter, ScalarFormatter, FormatStrFormatter, StrMethodFormatter
import matplotlib.pyplot as plt
import pandas as pd
import math
import numpy as np
import os
from enum import Enum
import re
from matplotlib.colors import LinearSegmentedColormap

BASELINE_TIME = .3
normalize_time = lambda n: n / BASELINE_TIME

class Mode(Enum):
    DEFAULT = "default"
    EVO = "evo"

class BenchmarkParams(Enum):
    ALGORITHM = "algorithm"
    GRAPH = "graph"
    SEED = "seed"
    K = "k"
    EPSILON = "epsilon"

class ResultParams(Enum):
    IMBALANCE = "imbalance"
    TIME = "totalPartitionTime"
    KM1 = "km1"
    CUT = "cut"
    FAILED = "failed"

class BenchmarkResult:
    def __init__(self, params: BenchmarkParams, results: ResultParams):
        self.params = params
        self.results = results

# Color map for difference matrix
DIFF_CMAP = LinearSegmentedColormap.from_list(
    "diff_red_green",
    [
        (0.00, "#300000"),
        (0.15, "#5a0000"),
        (0.30, "#b00000"),
        (0.45, "#ff4000"),
        (0.60, "#ff9800"),
        (0.75, "#e5ff00"),
        (0.90, "#6dff2b"),
        (1.00, "#00c400"),
    ],
    N=512
)


def parse_benchamark_input(input_file, mode: Mode):
    lowest_value = math.inf
    seed_index = 0
    with open(input_file, 'r', encoding='ascii') as file:
        lines = [line.strip() for line in file.readlines()]
        seed_count = lines.count('###') + 1
        seed_results = [[], []]
        max_results = 0
        results = 0
        all_seeds_results = [[] for _ in range(seed_count)]

        if (mode == Mode.EVO):
            for line in lines:
                if line == '###':
                    seed_index += 1
                    all_seeds_results[seed_index - 1] = seed_results
                    seed_results = [[], []]
                    results = 0
                    continue
                elif line.startswith('Starttime'):
                    start_time = datetime.fromtimestamp(int(line.split()[1])/1000)
                    continue
                timestamp, method, target_value = line.split(', ')
                timestamp = datetime.fromtimestamp(int(timestamp)/1000)
                target_value = int(target_value)
                delta = (timestamp - start_time).total_seconds()
                seed_results[0] += [normalize_time(delta)]
                seed_results[1] += [target_value]
                results += 1
                print(f"Time: {delta}, Method: {method}, Target Value: {target_value}")
            all_seeds_results[-1] = seed_results
            first_result = [(max(res[0][0] for res in all_seeds_results), np.mean([res[1][0] for res in all_seeds_results]))]
            result_array = first_result + [(math.inf,0) for i in range(sum(len(j[0])for j in all_seeds_results)-seed_count)] 
            print('Done')

        elif (mode == Mode.DEFAULT):
            start_time = datetime.strptime(lines[1][:-3], "%H:%M:%S.%f")
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
                    print('-'*20)
                    print(f"Start Time: {start_time} run {seed_index}")
                    continue
                timestamp, method, traget_value = line.split(', ')
                traget_value = int(traget_value)
                time = datetime.strptime(timestamp[:-3], "%H:%M:%S.%f")
                delta = None
                delta = (time - start_time).total_seconds()
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
    

def parse_runs_list_csv(input_file: str):
    """
    Parses a CSV file containing multiple runs.
    """



def parse_diff_matrices(diff_file: str, seed_id: int = 0):
    """
    Returns an array of difference matrices (each a 2D list) for a given seed_id.
    Matrices are separated by lines starting with ---.
    Seeds are separated by lines starting with ###.
    """
    seed_counter = 0
    matrices = []
    current = []
    with open(diff_file, 'r', encoding='ascii') as file:
        for raw in file:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('###'):
                if seed_counter == seed_id:
                    break 
                seed_counter += 1
                continue
            if seed_counter != seed_id:
                continue
            if line.startswith('---'):
                if current:
                    matrices.append(current)
                    current = []
                continue
            # row
            row = [float(x) for x in re.split(r',\s*', line) if x]
            current.append(row)
    if current:  # flush last
        matrices.append(current)
    return matrices


def plot_diff_matrix(diff_file: str,
                     seed_id: int = 0,
                     matrix_index: int | None = None,
                     cmap = DIFF_CMAP,
                     output_dir: str | None = None):
    """
    1) Saves a grid image containing all matrices for the given seed.
    2) Saves a single image for the matrix at matrix_index (if provided).
    3) Saves a line plot of average (off-diagonal) difference over time.
    """
    if output_dir is None:
        output_dir = os.path.dirname(os.path.abspath(__file__))

    matrices_raw = parse_diff_matrices(diff_file, seed_id=seed_id)
    if not matrices_raw:
        print(f"No diff matrices found for seed {seed_id} in {diff_file}")
        return

    matrices = [np.array(m) for m in matrices_raw]
    num = len(matrices)

    # 1) Grid of all matrices using each matrix' local max for normalization
    cols = math.ceil(math.sqrt(num))
    rows = math.ceil(num / cols)
    plt.figure(figsize=(3 * cols, 3 * rows))
    for i, m in enumerate(matrices):
        ax = plt.subplot(rows, cols, i + 1)
        local_max = int(m.max()) if m.size else 0
        if local_max == 0:
            local_max = 1
        im = ax.imshow(m, cmap=cmap, vmin=0, vmax=local_max, interpolation='nearest')
        ax.set_title(f"M{i}", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    plt.suptitle(f"Seed {seed_id} Diff Matrices (local normalization)", y=0.995)
    plt.tight_layout(rect=[0,0,1,0.97])
    all_path = os.path.join(output_dir, f"diff_seed{seed_id}_all.png")
    plt.savefig(all_path, dpi=160)
    plt.close()
    print(f"Saved all matrices grid -> {all_path}")

    # 2) Single matrix (local normalization so its own max is green)
    if matrix_index is not None:
        if 0 <= matrix_index < num:
            m = matrices[matrix_index]
            local_max = m.max()
            if local_max == 0:
                local_max = 1
            plt.figure(figsize=(5, 4))
            plt.imshow(m, cmap=cmap, vmin=0, vmax=local_max, interpolation='nearest')
            plt.title(f"Seed {seed_id} Matrix {matrix_index} (Max={int(local_max)})")
            plt.colorbar()
            plt.tight_layout()
            single_path = os.path.join(output_dir, f"diff_seed{seed_id}_matrix{matrix_index}.png")
            plt.savefig(single_path, dpi=160)
            plt.close()
            print(f"Saved matrix {matrix_index} -> {single_path}")
        else:
            print(f"matrix_index {matrix_index} out of range (0..{num-1})")

    # 3) Average off-diagonal over time
    avgs = []
    for m in matrices:
        if m.size == 0:
            avgs.append(0)
            continue
        mask = ~np.eye(m.shape[0], dtype=bool)
        off = m[mask]
        avgs.append(off.mean() if off.size else 0)
    plt.figure(figsize=(7, 4))
    plt.plot(range(num), avgs, marker='o', color="#b00000")
    plt.xlabel("Matrix index (time step)")
    plt.ylabel("Avg off-diagonal diff")
    plt.title(f"Average Difference Over Time (Seed {seed_id})")
    plt.grid(True, ls="--", alpha=0.6)
    avg_path = os.path.join(output_dir, f"diff_seed{seed_id}_avg_over_time.png")
    plt.tight_layout()
    plt.savefig(avg_path, dpi=160)
    plt.close()
    print(f"Saved avg diff plot -> {avg_path}")
    

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
    
    # plt.xscale('log')
    # plt.yscale('log')
    

    ax = plt.gca()

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.savefig(output_file)
    plt.close()

if __name__ == "__main__":
    # Multiple datasets:
    script_dir = os.path.dirname(os.path.abspath(__file__))

    result1 = parse_benchamark_input('output.csv', mode=Mode.DEFAULT)
    result2 = parse_benchamark_input('output_evo.csv', mode=Mode.EVO)
    datasets = [(result1, "Baseline"), (result2, "Evolutionary")]

    # Save output file in the same directory as the script
    output_path = os.path.join(script_dir, "comparison.png")
    create_plot(datasets, "Benchmark Comparison", "Time (t_n)", "Value", output_path, scientific=False)

    print(f"Plot saved to: {output_path}")

    diff_file_path = os.path.join(script_dir, "..", "output_evo_diff.csv")
    if os.path.exists(diff_file_path):
        plot_diff_matrix(diff_file_path, seed_id=0, matrix_index=0)
    else:
        print("Could not find diff file at", diff_file_path)
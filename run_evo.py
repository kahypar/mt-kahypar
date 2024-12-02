#!/usr/bin/env python

import argparse
import shlex
import os
import sys
import signal
import subprocess
from threading import Timer
from pathlib import Path

###########################


def get_args():
  parser = argparse.ArgumentParser()
  parser.add_argument("graph", type=str)
  parser.add_argument("--threads", type=int, required=True)
  parser.add_argument("--epsilon", type=float, default=0.03)
  parser.add_argument("--timelimit", type=int, required=True)  # for one k
  parser.add_argument("--maxk", type=int, default=32)

  return parser.parse_args()


def first_line(file):
    with open(file) as f:
        return f.readline().rstrip()


def run_mtk_evo(graph, timelimit, k, epsilon, threads, mt_kahypar, config, detect_instance_type=True, attempts=3):
  def kill_proc():
    os.killpg(os.getpgid(mt_kahypar_proc.pid), signal.SIGTERM)

  repetitions = 3
  path = Path(graph)
  freq_file = f"{path.stem}.freqk{k}.csv"
  freq_file = os.path.dirname(graph) + "/" + freq_file

  for i in range(0, attempts):
    if Path(freq_file).is_file():
      print(f"Skipping: {freq_file} {i}")
      return

    # Run MT-KaHyPar
    cmd = [mt_kahypar,
          "-h" + graph,
          "-p" + config,
          "-k" + str(k),
          "-e" + str(epsilon),
          "--seed=0",
          "-mdirect",
          "--s-num-threads=" + str(threads),
          "--verbose=false",
          "--evo-dynamic-time-limit=1",
          "--evo-randomized-flows=1",
          "--time-limit=" + str(int(timelimit / repetitions)),
          "--evo-repetitions=" + str(repetitions),
          "--evo-frequency-file=" + freq_file,
          ]
    if detect_instance_type:
      if graph.endswith(".metis") or graph.endswith(".graph"):
        cmd.append("--instance-type=graph")
        cmd.append("--input-file-format=metis")
        cmd.append("-ocut")
      else:
        cmd.append("-okm1")
    print(shlex.join(cmd))

    mt_kahypar_proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True, preexec_fn=os.setsid)

    t = Timer(timelimit + 1000, kill_proc)
    t.start()
    out, err = mt_kahypar_proc.communicate()
    t.cancel()

    print(out)
    if mt_kahypar_proc.returncode != 0:
      if err != None:
        print(err, file=sys.stderr)
    else:
      print(f"Success at attempt {i}")
      return
  print(f"Error: {freq_file}")


###################
# The main script #
###################

MT_KAHYPAR = "/home/nikolai/Documents/Hypergraphs/evo-mt-kahypar/build/mt-kahypar/application/MtKaHyPar"
CONFIG = "/home/nikolai/Documents/Hypergraphs/evo-mt-kahypar/config/evo_quality_preset.ini"

ks = [4, 8, 12, 16, 20, 24, 28, 32, 48, 64, 96, 128]

args = get_args()
for k in ks:
   timefactor = 1 + (k - 2) * 0.2
   timelimit = int(timefactor * args.timelimit)
   run_mtk_evo(args.graph, timelimit, k, args.epsilon, args.threads, MT_KAHYPAR, CONFIG)
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


def run_mtk_evo(graph, timelimit, k, epsilon, threads, mt_kahypar, config, detect_instance_type=True):
  repetitions = 5
  path = Path(graph)
  freq_file = f"{path.stem}.freqk{k}.csv"
  # Run MT-KaHyPar
  cmd = [mt_kahypar,
         "-h" + graph,
         "-p" + config,
         "-k" + str(k),
         "-e" + str(epsilon),
         "--seed=0",
         "-okm1",
         "-mdirect",
         "--s-num-threads=" + str(threads),
         "--verbose=false",
         "--time-limit=" + str(int(timelimit / repetitions)),
         "--evo-repetitions=" + str(repetitions),
         "--evo-frequency-file=" + os.path.dirname(graph) + "/" + freq_file,
        ]
  print(shlex.join(cmd))
  if detect_instance_type:
    if graph.endswith(".metis") or graph.endswith(".graph"):
      cmd.append("--instance-type=graph")
      cmd.append("--input-file-format=metis")

  mt_kahypar_proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True, preexec_fn=os.setsid)

  def kill_proc():
    os.killpg(os.getpgid(mt_kahypar_proc.pid), signal.SIGTERM)

  t = Timer(timelimit + 20, kill_proc)
  t.start()
  out, err = mt_kahypar_proc.communicate()
  t.cancel()

  if mt_kahypar_proc.returncode != 0:
    print(out)
    print(err, file=sys.stderr)



###################
# The main script #
###################

MT_KAHYPAR = "/home/nikolai/Documents/Hypergraphs/evo-mt-kahypar/build/mt-kahypar/application/MtKaHyPar"
CONFIG = "/home/nikolai/Documents/Hypergraphs/evo-mt-kahypar/config/evo_quality_preset.ini"

args = get_args()
for k in range(2, args.maxk + 1):
   # TODO: dependent on size of graph?
   timefactor = 1 + (k - 2) * 0.1
   timelimit = int(timefactor * args.timelimit)
   run_mtk_evo(args.graph, timelimit, k, args.epsilon, args.threads, MT_KAHYPAR, CONFIG)

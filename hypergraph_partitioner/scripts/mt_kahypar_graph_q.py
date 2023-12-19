#!/usr/bin/python3
import os
import os.path

from mt_kahypar_common import (get_args, invalid, parse, parse_or_default,
                               parse_required_value, print_result,
                               run_mtkahypar, set_result_vals)

###################################
# SETUP ENV
###################################
algorithm = "MT-KaHyPar-Graph-Q"
mt_kahypar = os.environ.get("MT_KAHYPAR")
assert (mt_kahypar != None), "check env.sh"
###################################

args = get_args()
if args.name != "":
  algorithm = args.name

result, success = run_mtkahypar(mt_kahypar, args, default_args={
  "--preset-type": "highest_quality",
  "--instance-type": "graph",
  "--input-file-format": "metis",
})

set_result_vals(
  km1=invalid,
  cut=invalid,
  total_time=invalid,
  imbalance=1.0,
)
if success:
  parse_required_value(result, "km1", parser=int)
  parse_required_value(result, "cut", parser=int)
  parse_required_value(result, "totalPartitionTime", out="total_time")
  parse_required_value(result, "imbalance")

# CSV format: algorithm,graph,timeout,seed,k,epsilon,num_threads,imbalance,totalPartitionTime,objective,km1,cut,failed
print_result(algorithm, args)

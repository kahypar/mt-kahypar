#!/bin/bash

EXEC="./build/mt-kahypar/application/MtKaHyPar.exe"
INSTANCE="./tests/instances/ibm01.hgr"

COMMON_ARGS="-h ${INSTANCE} -k 8 -e 0.03 -t 10 -m direct -okm1 --seed=${2:-1}"

extract_km1() {
  # extracts the km1 value
  grep -E '^.*km1\s+= [0-9]+ \(primary objective function\)' | awk '{print $3;}'
}

multilevel_run() {
  
  while true; do
    # Normal Multilevel (no evo)
    BASELINE_OUT="$(${EXEC} ${COMMON_ARGS} --preset-type=default 2>/dev/null | extract_km1)"
    TIMESTAMP=$(date +"%T.%N")
    echo "$TIMESTAMP, baseline, $BASELINE_OUT" >> output.csv
  done
}

evolutionary_run() {
  # Evolutionary enabled
  EVO_OUT="$(${EXEC} ${COMMON_ARGS} --preset-type=default --partition-evolutionary=true --time-limit=$1 2>/dev/null | extract_km1)"
  TIMESTAMP=$(date +"%T.%N")
  echo "$TIMESTAMP, evolutionary, $EVO_OUT" >> output_evo.csv
}

## MAIN
if [ -z "$1" ]; then
  echo "Usage: $0 mode"
  echo "mode: multilevel | evolutionary"
  exit 1
fi

if [ "$1" == "multilevel" ]; then
  multilevel_run 
elif [ "$1" == "evo" ]; then
  evolutionary_run "$3"
else
  echo "Invalid mode: $1"
  exit 1
fi
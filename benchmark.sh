#!/bin/bash

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]; then
  echo "Usage: $0 time_limit [seconds] seed_num [number of seeds] mode [evo || multilevel] history_file (required) [evo_diff_file (optional, only for evo mode)]"
  exit 1
fi

# help command
if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  echo "Usage: $0 time_limit [seconds] seed_num [number of seeds] mode [evo || multilevel] history_file (required) [evo_diff_file (optional, only for evo mode)]"
  exit 0
fi

if [ $# -lt 4 ]; then
  echo "Usage: $0 time_limit [seconds] seed_num [number of seeds] mode [evo || multilevel] history_file (required) [evo_diff_file (optional, only for evo mode)]"
  exit 1
fi

EXEC="./build/mt-kahypar/application/MtKaHyPar.exe"
INSTANCE="./tests/instances/ibm01.hgr"

COMMON_ARGS="-h ${INSTANCE} -k 8 -e 0.03 -t 10 -m direct -okm1"

TIME_LIMIT="$1"
SEED_NUM="$2"
MODE="$3"
# 4th arg: history output file (required for both modes)
# 5th arg: diff matrix output file (optional, only for evo mode)
HISTORY_FILE="$4"
EVO_DIFF_FILE="$5"



if [ "$MODE" == "multilevel" ]; then
  echo "timestamp, mode, km1" > "${HISTORY_FILE}"

  # loop for every seed
  for i in $(seq 1 "$SEED_NUM"); do
  START_TIMESTAMP=$(date +"%T.%N")
  echo "$START_TIMESTAMP" >> "${HISTORY_FILE}"
    # echo "Running iteration $i with seed $i (${TIME_LIMIT}s timeout)"
    timeout "$TIME_LIMIT" ./infinite-run.sh "$MODE" "$i"
    if [ "$i" -lt "$SEED_NUM" ]; then
      echo "###" >> "${HISTORY_FILE}"
    fi
  done

elif [ "$MODE" == "evo" ]; then

  # Ensure history file (required)
  if [ -z "$HISTORY_FILE" ]; then
    echo "History output file is required as the 4th argument."
    echo "Usage: $0 time_limit seed_num evo|multilevel history_file [evo_diff_file]"
    exit 1
  fi

  # prepare flags
  EVO_DIFF_FLAG=""
  if [ -n "$EVO_DIFF_FILE" ]; then
    EVO_DIFF_FLAG="--evo-diff-matrix-file=${EVO_DIFF_FILE}"
    echo "Evo diff matrix file set to: $EVO_DIFF_FILE"
    # clear evo diff file
    > "$EVO_DIFF_FILE"
  fi

  # clear history file
  > "$HISTORY_FILE"
  echo "History file set to: $HISTORY_FILE"

  # loop for every seed
  for i in $(seq 1 "$SEED_NUM"); do
    echo "Running iteration $i with seed $i (${TIME_LIMIT}s timeout)"

    ${EXEC} ${COMMON_ARGS} --preset-type=default \
      --seed=$i --partition-evolutionary=true \
      --time-limit=${TIME_LIMIT} --enable-benchmark=1 \
      --evo-history-file="${HISTORY_FILE}" ${EVO_DIFF_FLAG} 2>/dev/null

    if [ "$i" -lt "$SEED_NUM" ]; then
      echo "###" >> "${HISTORY_FILE}"
      if [ -n "$EVO_DIFF_FILE" ]; then
        echo "###" >> "${EVO_DIFF_FILE}"
      fi
    fi
  done

else
  echo "Invalid mode: $MODE"
  exit 1
fi

extract_km1() {
  # extracts the km1 value
  grep -E '^.*km1\s+= [0-9]+ \(primary objective function\)' | awk '{print $3;}'
}

multilevel_run() {
  
  while true; do
    # Normal Multilevel (no evo)
    BASELINE_OUT="$(${EXEC} ${COMMON_ARGS} --preset-type=default 2>/dev/null | extract_km1)"
      TIMESTAMP=$(date +"%T.%N")
      echo "$TIMESTAMP, baseline, $BASELINE_OUT" >> "${HISTORY_FILE}"
  done
}

evolutionary_run() {
  EVO_OUT="$(${EXEC} ${COMMON_ARGS} --preset-type=default --partition-evolutionary=true --time-limit=$1 2>/dev/null | extract_km1)"
  echo "$TIMESTAMP, evolutionary, $EVO_OUT" >> output_evo.csv
}
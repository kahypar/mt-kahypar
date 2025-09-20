#!/bin/bash

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]; then
  echo "Usage: $0 time_limit [seconds] seed_num [number of seeds] mode [evo || multilevel] time_step [time step limit for evo]"
  exit 1
fi

# help command
if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  echo "Usage: $0 time_limit [seconds] seed_num [number of seeds] mode [evo || multilevel]"
  exit 0
fi

TIME_LIMIT="$1"
SEED_NUM="$2"
MODE="$3"
TIME_STEP="$4"


if [ "$MODE" == "multilevel" ]; then
  echo "timestamp, mode, km1" > output.csv

  # loop for every seed
  for i in $(seq 1 "$SEED_NUM"); do
    START_TIMESTAMP=$(date +"%T.%N")
    echo "$START_TIMESTAMP" >> output.csv
    # echo "Running iteration $i with seed $i (${TIME_LIMIT}s timeout)"
    timeout "$TIME_LIMIT" ./infinite-run.sh "$MODE" "$i"
    if [ "$i" -lt "$SEED_NUM" ]; then
      echo "###" >> output.csv
    fi
  done
elif [ "$MODE" == "evo" ]; then
  echo "timestamp, mode, km1" > output_evo.csv

  # loop for every seed
  for i in $(seq 1 "$SEED_NUM"); do
    echo "Running seed $i"
    START_TIMESTAMP=$(date +"%T.%N")
    echo "$START_TIMESTAMP" >> output_evo.csv
    # loop for every time step multiple
    for j in $(seq 1 "$((TIME_LIMIT / TIME_STEP))"); do
      current_time_limit=$((j * TIME_STEP))
      echo "Running iteration $i with seed $i (${current_time_limit}s timeout)"
      ./infinite-run.sh "$MODE" "$i" "$current_time_limit"
    done
    if [ "$i" -lt "$SEED_NUM" ]; then
      echo "###" >> output_evo.csv
    fi
  done

else
  echo "Invalid mode: $MODE"
  exit 1
fi
#!/bin/bash

# Check if argument is provided
if [ -z "$1" ] || [ -z "$2" ]; then
  echo "Usage: $0 n [runs] t [seconds per run]"
  exit 1
fi

# help command
if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  echo "Usage: $0 n [runs] t [seconds per run]"
  exit 0
fi

echo "timestamp, km1" > output.csv


for i in $(seq 1 "$1"); do
  echo "Second: $i"
  timeout "$2" ./infinite-run.sh
  if [ "$i" -lt "$1" ]; then
    echo "###" >> output.csv
  fi
done

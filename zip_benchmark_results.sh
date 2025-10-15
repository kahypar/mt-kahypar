#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: $0 <benchmark_directory>"
    exit 1
fi

BENCHMARK_DIR="$1"

# Check if benchmark directory exists
if [ ! -d "$BENCHMARK_DIR" ]; then
    echo "Error: Benchmark directory '$BENCHMARK_DIR' not found!"
    exit 1
fi

BENCHMARK_NAME=$(basename "$BENCHMARK_DIR")

ZIP_FILE="${BENCHMARK_DIR}/${BENCHMARK_NAME}_results.zip"

echo "Zipping results from: $BENCHMARK_DIR"

cd "$BENCHMARK_DIR" || exit 1
zip -r "${BENCHMARK_NAME}_results.zip" \
    *.csv \
    *_results/ \
    workload.txt \
    -x "*.header.csv"

echo "Results zipped to: ${ZIP_FILE}"

# Print zip contents
echo ""
echo "Archive contents:"
unzip -l "${BENCHMARK_NAME}_results.zip"
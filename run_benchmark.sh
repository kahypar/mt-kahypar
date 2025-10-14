#!/bin/bash
// filepath: run_benchmark.sh

# Check if config file is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <config.json>"
    echo "Example: $0 benchmarks/config.json"
    exit 1
fi

CONFIG_FILE="$1"

# Check if config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file '$CONFIG_FILE' not found!"
    exit 1
fi

echo "Setting up benchmark with config: $CONFIG_FILE"

python3 benchmarks/setup_benchmark.py "$CONFIG_FILE" -f

if [ $? -ne 0 ]; then
    echo "Error: Benchmark setup failed!"
    exit 1
fi

echo "Setup completed successfully. Starting benchmark execution..."

python3 benchmarks/execute_benchmark.py "$CONFIG_FILE"

if [ $? -ne 0 ]; then
    echo "Error: Benchmark execution failed!"
    exit 1
fi

echo "Benchmark completed successfully!"
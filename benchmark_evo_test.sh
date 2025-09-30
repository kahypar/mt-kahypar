#!/bin/bash

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
  echo "Usage: $0 k [different runs] n [seconds per run] output_file"
  exit 1
fi

extract_km1() {
  grep -E '^.*km1\s+= [0-9]+ \(primary objective function\)' | awk '{print $3;}'
}

EXEC="./build/mt-kahypar/application/MtKaHyPar.exe"
INSTANCE="./tests/instances/ibm01.hgr"

COMMON_ARGS="-h ${INSTANCE} -k 8 -e 0.03 -t 10 -m direct --objective=km1"

K="$1"
N="$2"
OUTPUT_FILE="$3"

if [ -n "$OUTPUT_FILE" ]; then
  echo "Output file set to: $OUTPUT_FILE"
  # clear output file
  > "$OUTPUT_FILE"
fi

# multiple short runs
echo "multiple short runs results: " >> "${OUTPUT_FILE}"
for i in $(seq 1 "$K"); do
    echo "Running iteration $i with seed $i (${N}s timeout)"

    output_short=$(${EXEC} ${COMMON_ARGS} --preset-type=default \
        --seed=$i --partition-evolutionary=true \
        --time-limit=${N} 2>/dev/null | extract_km1)
    echo "output short: ${output_short}"
    echo "run $i: ${output_short}" >> "${OUTPUT_FILE}"
done

# singular long run (N*K seconds)
LONG_RUN_SECONDS=$((N * K))
echo "-----" >> "${OUTPUT_FILE}" 
output_long=$(${EXEC} ${COMMON_ARGS} --preset-type=default \
        --seed=1 --partition-evolutionary=true \
        --time-limit=${LONG_RUN_SECONDS} 2>/dev/null | extract_km1)

echo "long run: ${output_long}" >> "${OUTPUT_FILE}"
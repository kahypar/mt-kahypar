#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: $0 <benchmark_directory_or_prefix>"
    exit 1
fi

INPUT="$1"

# Function to zip a single benchmark directory
zip_single_benchmark() {
    local BENCHMARK_DIR="$1"
    local BENCHMARK_NAME=$(basename "$BENCHMARK_DIR")
    local ZIP_FILE="${BENCHMARK_DIR}/${BENCHMARK_NAME}_results.zip"

    echo "Zipping results from: $BENCHMARK_DIR"

    cd "$BENCHMARK_DIR" || return 1
    zip -r "${BENCHMARK_NAME}_results.zip" \
        *.csv \
        *_results/ \
        workload.txt \
        -x "*.header.csv"

    echo "Results zipped to: ${ZIP_FILE}"
    echo ""
    echo "Archive contents:"
    unzip -l "${BENCHMARK_NAME}_results.zip"
    echo ""
    
    cd - > /dev/null
}

# Check if input is an exact directory
if [ -d "$INPUT" ]; then
    zip_single_benchmark "$INPUT"
else
    # Find all matching directories
    # Extract the directory path and the prefix
    DIR_PATH=$(dirname "$INPUT")
    PREFIX=$(basename "$INPUT")
    
    # Search for directories matching the prefix in the specified directory
    MATCHING_DIRS=$(find "$DIR_PATH" -maxdepth 1 -type d -name "${PREFIX}*" 2>/dev/null | sort)
    
    if [ -z "$MATCHING_DIRS" ]; then
        echo "Error: No directories found matching prefix '$PREFIX' in '$DIR_PATH'"
        exit 1
    fi
    
    # Count matching directories
    DIR_COUNT=$(echo "$MATCHING_DIRS" | wc -l)
    
    echo "Found $DIR_COUNT director(ies) matching prefix '$PREFIX' in '$DIR_PATH':"
    echo "$MATCHING_DIRS"
    echo ""
    
    # Create combined zip file
    COMBINED_ZIP="${DIR_PATH}/${PREFIX}_combined_results.zip"
    
    # Build array of directory names (basenames only)
    DIR_NAMES=()
    while IFS= read -r dir; do
        if [ -d "$dir" ]; then
            DIR_NAMES+=("$(basename "$dir")")
        fi
    done <<< "$MATCHING_DIRS"
    
    # Zip all matching directories together
    cd "$DIR_PATH" || exit 1
    
    # Create the zip file with all directories
    zip -r "${PREFIX}_combined_results.zip" "${DIR_NAMES[@]}"
    
    cd - > /dev/null
    
    echo ""
    echo "Results zipped to: ${COMBINED_ZIP}"
    echo ""
    echo "Archive contents:"
    unzip -l "${COMBINED_ZIP}"
    echo ""
    echo "Summary: Zipped $DIR_COUNT benchmark directories into one archive"
fi
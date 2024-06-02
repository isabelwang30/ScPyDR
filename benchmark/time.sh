#!/bin/bash

# Make this script executable
chmod +x "$0"  # "$0" represents the path to the current script

# Define input variables
DATADIR="benchmark"
OUTPUT_DIR="scpydr_results"

# Run scpydr PCA
echo "Running scpydr PCA..."
mkdir -p $OUTPUT_DIR
time scpydr benchmark -o $OUTPUT_DIR

# Run scanpy PCA
echo "Running scanpy PCA..."
time python scanpy_pca_script.py $DATADIR -o $OUTPUT_DIR

# Compare results if needed
# Add your comparison logic here
#!/bin/bash

# Make this script executable
chmod +x "$0"  # "$0" represents the path to the current script

# Define input variables
DATADIR="benchmark"

# Run scpydr PCA
echo "Running scpydr PCA..."
# get time output while discarding scpydr terminal output
echo "Total scpydr PCA time:"
time scpydr $DATADIR -o $DATADIR/scpydr_results > /dev/null 2>&1
echo  # Insert empty line after scpydr time output

# Run scanpy PCA
echo "Running scanpy PCA..."
echo "Total scanpy PCA time:"
time python $DATADIR/scanpy_pca.py $DATADIR -o $DATADIR/scanpy_results
# Compare results if needed
# Add your comparison logic here

#!/bin/bash
# run with bash benchmark/time.sh

# Make this script executable
chmod +x "$0"  # "$0" represents the path to the current script

# Define input variables
BENCHDIR="benchmark"
DATADIR="benchmark/data"

# -------------------- PCA --------------------
# Run scpydr PCA
echo "Running scpydr PCA..."
# get time output while discarding scpydr terminal output
echo "Total scpydr PCA time:"
# run with 10 PCs for comparison
time scpydr $DATADIR -o $BENCHDIR/scpydr_results -n 10 > /dev/null 2>&1
echo  # Insert empty line after time output

# Run scanpy PCA
echo "Running scanpy PCA..."
echo "Total scanpy PCA time:"
# running with 10 PCs for comparison
time python $BENCHDIR/scanpy_pca.py $DATADIR -o $BENCHDIR/scanpy_results > /dev/null 2>&1
echo  # Insert empty line after time output

# -------------------- UMAP --------------------
# Run scpydr UMAP
echo "Running scpydr UMAP..."
echo "Total scpydr UMAP time:"
time python $BENCHDIR/scpydr_umap.py $DATADIR -o $BENCHDIR/scpydr_results > /dev/null 2>&1
echo  # Insert empty line after time output

# Run scanpy UMAP
echo "Running scanpy UMAP..."
echo "Total scanpy UMAP time:"
time python $BENCHDIR/scanpy_umap.py $DATADIR -o $BENCHDIR/scanpy_results > /dev/null 2>&1
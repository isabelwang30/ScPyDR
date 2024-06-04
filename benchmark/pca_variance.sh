#!/bin/bash
# run with bash benchmark/pca_variance.sh

# Make this script executable
chmod +x "$0"  # "$0" represents the path to the current script

# Define input variables
DATADIR="benchmark/data"


# Run scpydr PCA
echo "Running scpydr PCA..."
python benchmark/scpydr_pca_variance.py $DATADIR -n 10
echo  # Insert empty line after scpydr output


# Run scanpy PCA
echo "Running scanpy PCA..."
python benchmark/scanpy_pca_variance.py $DATADIR -n 10

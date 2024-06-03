#!/bin/bash


# Make this script executable
chmod +x "$0"  # "$0" represents the path to the current script

# Define input variables
DATADIR="data"


# Run scpydr PCA
echo "Running scpydr PCA..."
python3 scpydr_pca_variance.py $DATADIR -n 10
echo  # Insert empty line after scpydr output


# Run scanpy PCA
echo "Running scanpy PCA..."
python3 scanpy_pca_variance.py $DATADIR -n 10
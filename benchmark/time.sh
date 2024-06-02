# run with command: bash benchmark/time.sh
#!/bin/bash

# Make this script executable
chmod +x "$0"  # "$0" represents the path to the current script

# Define input variables
DATADIR="benchmark"

# Run scpydr PCA
echo "Running scpydr PCA..."
TOTAL_SCPYDR_TIME=$(echo $( TIMEFORMAT="%3U + %3S"; { time scpydr $DATADIR -o $DATADIR/scpydr_results; } 2>&1) | bc -l)
echo "Total scpydr PCA time: $TOTAL_SCPYDR_TIME seconds"

# Run scanpy PCA
echo "Running scanpy PCA..."
TOTAL_SCANPY_TIME=$(echo $( TIMEFORMAT="%3U + %3S"; { time python $DATADIR/scanpy_pca.py $DATADIR -o $DATADIR/scanpy_results; } 2>&1) | bc -l)
echo "Total scanpy PCA time: $TOTAL_SCANPY_TIME seconds"

# Compare results if needed
# Add your comparison logic here

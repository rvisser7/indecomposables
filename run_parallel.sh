#!/usr/bin/env bash

set -e

DEGREE=$1
DISC_MIN=$2
DISC_MAX=$3
VERBOSE=$4
NUM_THREADS=$5

if [ -z "$DEGREE" ] || [ -z "$DISC_MAX" ] || [ -z "$NUM_THREADS" ]; then
    echo "Usage: $0 DEGREE DISC_MIN DISC_MAX VERBOSE NUM_THREADS"
    echo "Example: $0 3 1 10000 1 20"
    exit 1
fi

echo "Launching parallel indecomposable computations"
echo "Degree: $DEGREE"
echo "Discriminant range: ${DISC_MIN:-1} to $DISC_MAX"
echo "Verbose: ${VERBOSE:-0}"
echo "Number of threads: $NUM_THREADS"

for ((i=0; i < NUM_THREADS; i=i+1))
do
    session_name="indecomposables_${DEGREE}_${DISC_MIN:-1}_${DISC_MAX}_${i}"
    cmd="./compute_indecomposables.sh $DEGREE ${DISC_MIN:-1} $DISC_MAX ${VERBOSE:-0} $NUM_THREADS $i"
    screen -S "$session_name" -dm bash -lc "$cmd"
    echo "Started screen session $session_name"
done

echo "All sessions launched. Use 'screen -ls' to view them."

#!/usr/bin/env bash

set -e

DEGREE=$1
DISC_MIN=$2
DISC_MAX=$3
VERBOSE=$4
NUM_THREADS=${5:-1}
THREAD_ID=${6:-0}
DEBUG=${7:-0}

if [ -z "$DEGREE" ] || [ -z "$DISC_MAX" ]; then
  echo "Usage: $0 DEGREE DISC_MIN DISC_MAX [VERBOSE] [NUM_THREADS] [THREAD_ID] [DEBUG]"
  echo "VERBOSE: optional, set to 1 to enable verbose output"
  echo "NUM_THREADS: optional, default 1"
  echo "THREAD_ID: optional, default 0"
  echo "DEBUG: optional, set to 1 to enable debug logs in console"
  exit 1
fi

echo "Running big indecomposable computation..."
echo "Degree: $DEGREE"
echo "Discriminant range: $DISC_MIN to $DISC_MAX"
echo "Threads: $NUM_THREADS"
echo "Thread ID: $THREAD_ID"
echo "Debug: ${DEBUG}"
if [ "$VERBOSE" = "1" ]; then
  echo "Verbose: enabled"
  VERBOSE_FLAG="--verbose"
else
  echo "Verbose: disabled"
  VERBOSE_FLAG=""
fi

if [ "$DEBUG" = "1" ]; then
  echo "Debug logging: enabled"
  DEBUG_FLAG="--debug"
else
  echo "Debug logging: disabled"
  DEBUG_FLAG=""
fi

sage -python scripts/compute_indecomposables.py \
  --degree $DEGREE \
  --disc-min ${DISC_MIN:-1} \
  --disc-max $DISC_MAX \
  --threads $NUM_THREADS \
  --thread-id $THREAD_ID \
  $DEBUG_FLAG \
  $VERBOSE_FLAG

echo "Done!"

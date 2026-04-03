#!/usr/bin/env bash

set -e

DEGREE=$1
DISC_MIN=$2
DISC_MAX=$3
VERBOSE=$4

if [ -z "$DEGREE" ] || [ -z "$DISC_MAX" ]; then
  echo "Usage: $0 DEGREE DISC_MIN DISC_MAX [VERBOSE]"
  echo "VERBOSE: optional, set to 1 to enable verbose output"
  exit 1
fi

echo "Running indecomposable computation..."
echo "Degree: $DEGREE"
echo "Discriminant range: $DISC_MIN to $DISC_MAX"
if [ "$VERBOSE" = "1" ]; then
  echo "Verbose: enabled"
  VERBOSE_FLAG="--verbose"
else
  echo "Verbose: disabled"
  VERBOSE_FLAG=""
fi

sage -python scripts/compute_indecomposables.py \
  --degree $DEGREE \
  --disc-min ${DISC_MIN:-1} \
  --disc-max $DISC_MAX \
  --output results_deg${DEGREE}.txt \
  $VERBOSE_FLAG

echo "Done!"

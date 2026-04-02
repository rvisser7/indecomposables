#!/usr/bin/env bash

set -e

DEGREE=$1
DISC_MIN=$2
DISC_MAX=$3

if [ -z "$DEGREE" ] || [ -z "$DISC_MAX" ]; then
  echo "Usage: $0 DEGREE DISC_MIN DISC_MAX"
  exit 1
fi

echo "Running indecomposable computation..."
echo "Degree: $DEGREE"
echo "Discriminant range: $DISC_MIN to $DISC_MAX"

sage -python scripts/compute_indecomposables.py \
  --degree $DEGREE \
  --disc-min ${DISC_MIN:-1} \
  --disc-max $DISC_MAX \
  --output results_deg${DEGREE}.txt

echo "Done!"

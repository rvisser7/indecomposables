#!/usr/bin/env bash

set -e

echo "Running real quadratic computation..."

sage -python scripts/compute_real_quadratics.py \
  --D-max 10000 \
  --output results_quadratic.txt

echo "Done!"


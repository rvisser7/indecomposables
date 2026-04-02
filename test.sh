#!/usr/bin/env bash
# This script runs pytest for the whole indecomposables project.

set -e  # stop on first error

echo "=============================="
echo " Running test suite"
echo "=============================="

echo "Running tests with Sage..."

echo "Running pytest..."
sage -python -m pytest tests/ -v

echo "=============================="
echo " All tests passed!"
echo "=============================="



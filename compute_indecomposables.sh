#!/usr/bin/env bash

set -e

print_usage() {
  echo "Usage: $0 DEGREE DISC_MIN DISC_MAX [options]"
  echo ""
  echo "Options:"
  echo "  -v, --verbose          Enable verbose console output"
  echo "  -d, --debug            Enable debug console logging"
  echo "  --num-threads N        Set total thread count (default: 1)"
  echo "  --thread-id I          Set thread id (default: 0)"
  echo "  -h, --help             Show this help"
}

if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
  print_usage
  exit 0
fi

DEGREE=${1:-}
DISC_MIN=${2:-}
DISC_MAX=${3:-}

if [[ $# -ge 3 ]]; then
  shift 3
else
  set --
fi

# Defaults for optional settings
VERBOSE=0
DEBUG=0
NUM_THREADS=1
THREAD_ID=0

if [ -z "$DEGREE" ] || [ -z "$DISC_MAX" ]; then
  print_usage
  exit 1
fi

# Named options
while [[ $# -gt 0 ]]; do
  case "$1" in
    -v|--verbose)
      VERBOSE=1
      shift
      ;;
    -d|--debug)
      DEBUG=1
      shift
      ;;
    --num-threads)
      NUM_THREADS="$2"
      shift 2
      ;;
    --thread-id)
      THREAD_ID="$2"
      shift 2
      ;;
    -h|--help)
      print_usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      print_usage
      exit 1
      ;;
  esac
done

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

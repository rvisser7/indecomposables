#!/usr/bin/env bash

set -e

print_usage() {
    echo "Usage: $0 DEGREE DISC_MIN DISC_MAX NUM_THREADS [options]"
    echo ""
    echo "Options:"
    echo "  -v, --verbose          Enable verbose mode for each worker"
    echo "  -h, --help             Show this help"
    echo ""
    echo "Example: $0 3 1 10000 20 -v"
}

if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
    print_usage
    exit 0
fi

DEGREE=${1:-}
DISC_MIN=${2:-}
DISC_MAX=${3:-}
NUM_THREADS=${4:-}

if [[ $# -ge 4 ]]; then
    shift 4
else
    set --
fi

VERBOSE=0

if [ -z "$DEGREE" ] || [ -z "$DISC_MAX" ] || [ -z "$NUM_THREADS" ]; then
    print_usage
    exit 1
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        -v|--verbose)
            VERBOSE=1
            shift
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

echo "Launching parallel indecomposable computations"
echo "Degree: $DEGREE"
echo "Discriminant range: ${DISC_MIN:-1} to $DISC_MAX"
echo "Verbose: $VERBOSE"
echo "Number of threads: $NUM_THREADS"

for ((i=0; i < NUM_THREADS; i=i+1))
do
    session_name="compute_indecomposables_${DEGREE}_${DISC_MIN:-1}_${DISC_MAX}_${NUM_THREADS}_${i}"
    cmd="./compute_indecomposables.sh $DEGREE ${DISC_MIN:-1} $DISC_MAX --num-threads $NUM_THREADS --thread-id $i"
    if [ "$VERBOSE" = "1" ]; then
        cmd="$cmd -v"
    fi
    screen -S "$session_name" -dm bash -lc "$cmd"
    echo "Started screen session $session_name"
done

echo "All sessions launched. Use 'screen -ls' to view them."

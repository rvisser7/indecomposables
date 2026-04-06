#!/usr/bin/env bash

set -e

print_usage() {
    echo "Usage: $0 [DEGREE DISC_MIN DISC_MAX NUM_THREADS]"
    echo ""
    echo "Behavior:"
    echo "  With no arguments: close all compute_indecomposables_* screen sessions"
    echo "  With 4 arguments: close only sessions from that specific run"
    echo ""
    echo "Examples:"
    echo "  $0"
    echo "  $0 3 1 10000 20"
}

if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
    print_usage
    exit 0
fi

# If args are provided, stop a specific run started by run_parallel.sh:
#   end_parallel.sh DEGREE DISC_MIN DISC_MAX NUM_THREADS
# If no args are provided, stop all sessions whose name starts with "compute_indecomposables_".

DEGREE=$1
DISC_MIN=$2
DISC_MAX=$3
NUM_THREADS=$4

close_session_if_exists() {
    session_name="$1"
    if screen -ls | grep -q "[.]${session_name}[[:space:]]"; then
        screen -S "$session_name" -X quit
        echo "Closed screen session $session_name"
    fi
}

if [ -z "$DEGREE" ] && [ -z "$DISC_MIN" ] && [ -z "$DISC_MAX" ] && [ -z "$NUM_THREADS" ]; then
    echo "Closing all compute_indecomposables screen sessions..."
    sessions=$(screen -ls | grep 'compute_indecomposables_' | awk '{print $1}' || true)

    if [ -z "$sessions" ]; then
        echo "No compute_indecomposables screen sessions found."
        exit 0
    fi

    while IFS= read -r full_id; do
        [ -z "$full_id" ] && continue
        session_name=${full_id#*.}
        close_session_if_exists "$session_name"
    done <<< "$sessions"

        echo "Done closing compute_indecomposables sessions."
    exit 0
fi

if [ -z "$DEGREE" ] || [ -z "$DISC_MAX" ] || [ -z "$NUM_THREADS" ]; then
        print_usage
    exit 1
fi

echo "Closing compute_indecomposables sessions for degree=$DEGREE, discs=${DISC_MIN:-1}..$DISC_MAX, threads=$NUM_THREADS"
for ((i=0; i < NUM_THREADS; i=i+1))
do
    session_name="compute_indecomposables_${DEGREE}_${DISC_MIN:-1}_${DISC_MAX}_${NUM_THREADS}_${i}"
    close_session_if_exists "$session_name"
done

echo "Done."

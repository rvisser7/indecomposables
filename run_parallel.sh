#!/bin/bash

# A simple bash script to run compute_indecomposables in parallel 
# across multiple screen sessions (similar to genus2 code)

# Running threads in parallel
for ((i=1; i <= 7; i=i+1))
do
    screen -S indecomposables_${i} -md nice compute_indecomposables.sh ${i}
done

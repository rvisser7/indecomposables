#!/bin/bash
# Shell script to automatically end all screen sessions running indecomposables code in parallel!

# Ending all threads in parallel
for ((i=1; i <= 7; i=i+1))
do
    # End the KL field system
    screen -S indecomposables_${i} -X quit
done

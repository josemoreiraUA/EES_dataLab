#!/bin/bash

# Loop through dist_range values 1 to 10
for dist in {1..10}
do
    # Call Python program with current dist_range value
    python tests.py $dist

    # Append dist_range and totalPoints to output file
    echo "$dist,$(tail -1 output.txt)" >> results.txt
done



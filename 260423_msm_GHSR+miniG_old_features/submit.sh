#!/bin/bash

base=$(pwd)

for seed in */; do
    seed=${seed%/}
    for j in {0..3}; do
        cd "$base/$seed/$j" || continue
        cp "$base/discard.sh" .
        pwd
        sbatch discard.sh
    done
done



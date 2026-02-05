#!/bin/bash

base=$(pwd)
declare -a systems=( human mouse )

for system in "${systems[@]}";do
    cd "$base/$system" || continue
    for binder in */; do
        binder=${binder%/}
        for j in {0..2}; do
            cd "$base/$system/$binder/$j" || continue
            cp "$base/job.sh" .
            pwd
            sbatch job.sh
        done
    done
done


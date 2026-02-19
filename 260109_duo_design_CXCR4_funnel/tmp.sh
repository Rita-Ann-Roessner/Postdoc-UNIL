#!/bin/bash

pwd=$(pwd)

declare -a params=( human mouse )

for batch in "${params[@]}"
do
    cd $batch
    for dir in */
    do
        binder=${dir%/}
        for i in {0..2}
        do
            dir="$pwd/$batch/$binder/$i"
            cd $dir
            plumed sum_hills --hills HILLS --outfile fes.dat 
            # every 10 ns
            rm -r FES
            mkdir FES
            plumed sum_hills --hills HILLS --outfile FES/fes. --stride 10000 --mintozero 
        done
    done
    cd $pwd
done

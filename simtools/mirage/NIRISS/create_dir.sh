#!/bin/bash

for folder in image378_sim image476_sim image430_sim image277_sim image378_data image476_data image430_data image277_data
do
    cd ${folder}
    for dir in 10jup50au 5jup50au 1jup50au 1sat50au 10jup30au 5jup30au 1jup30au
    do
        cd ${dir}
        mkdir 5
        mkdir 0
        mkdir 1
        mkdir 2
        cd ..
    done
    cd ..
done
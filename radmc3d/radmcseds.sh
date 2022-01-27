#!/bin/bash

sed_nostar(){
cd $1
cd $2
radmc3d sed sizeau $3 fluxcons nostar
cp spectrum.out sed_nostar.out
cd ..
cd ..
}


sed_star(){
cd $1
cd $2
radmc3d sed sizeau $3 fluxcons
cp spectrum.out sed_star.out
cd ..
cd ..
}

for dir in 10jup50au 5jup50au 1jup50au 1sat50au 10jup30au 5jup30au 1jup30au 
do
    for casedir in ori 0 1 2 3 4
    do
        sed_nostar ${dir} ${casedir} 250 &
    done
done

for dir in 10jup50au 5jup50au 1jup50au 1sat50au 10jup30au 5jup30au 1jup30au 
do
    for casedir in 2
    do
        sed_nostar ${dir} ${casedir} 250 &
    done
done

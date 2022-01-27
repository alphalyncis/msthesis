#!/bin/bash


# lam 378
for dir in 10jup50au 5jup50au 1jup50au 1sat50au 
do
    for casedir in 2
    do
        for inc in 0 30 60
        do
            python metis_sim.py ${dir}/${casedir}/npix520lam378inc${inc}.fits coron
        done
    done
done

for dir in 10jup30au 5jup30au 1jup30au
do
    for casedir in 2
    do
        for inc in 0 30 60
        do
            python metis_sim.py ${dir}/${casedir}/npix300lam378inc${inc}.fits coron
        done
    done
done

# lam 476
for dir in 10jup50au 5jup50au 1jup50au 1sat50au 
do
    for casedir in 2
    do
        for inc in 0 30 60
        do
            python metis_sim.py ${dir}/${casedir}/npix520lam476inc${inc}.fits coron
        done
    done
done

for dir in 10jup30au 5jup30au 1jup30au
do
    for casedir in 2
    do
        for inc in 0 30 60
        do
            python metis_sim.py ${dir}/${casedir}/npix300lam476inc${inc}.fits coron
        done
    done
done


#lam 1020
for dir in 10jup50au 5jup50au 1jup50au 1sat50au 
do
    for casedir in 2
    do
        for inc in 0 30 60
        do
            python metis_sim.py ${dir}/${casedir}/npix371lam1020inc${inc}.fits coron
        done
    done
done

for dir in  1jup30au 10jup30au 5jup30au
do
    for casedir in 2
    do
        for inc in 0 30 60
        do
            python metis_sim.py ${dir}/${casedir}/npix214lam1020inc${inc}.fits coron
        done
    done
done

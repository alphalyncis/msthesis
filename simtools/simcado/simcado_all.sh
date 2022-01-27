#!/bin/bash
for lam in 214 125 164
do

    for dir in 10jup50au 5jup50au 1jup50au 1sat50au 
    do
        for casedir in 2
        do
            for inc in 0 30 60
            do
                python micado_sim.py ${dir}/${casedir}/npix1733lam${lam}inc${inc}.fits coron
            done
        done
    done

    for dir in 10jup30au 5jup30au  1jup30au
    do
        for casedir in 2
        do
            for inc in 0 30 60
            do
                python micado_sim.py ${dir}/${casedir}/npix1000lam${lam}inc${inc}.fits coron
            done
        done
    done

done
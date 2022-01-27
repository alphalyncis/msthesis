#!/bin/bash

for lam in 1020 1500 2100
do
    for dir in 10jup50au 5jup50au 1jup50au 1sat50au 
    do
        for casedir in 0 1 2
        do
            for inc in 0 30 60
            do
                python coron_mask.py ${dir}/${casedir}/cube5npix23lam${lam}inc${inc}.fits
            done
        done
    done

    for dir in 10jup30au 5jup30au 1jup30au
    do
        for casedir in 0 1 2
        do
            for inc in 0 30 60
            do
                python coron_mask.py ${dir}/${casedir}/cube5npix13lam${lam}inc${inc}.fits
            done
        done
    done
done
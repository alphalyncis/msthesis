#!/bin/bash

declare -a shortfilters=(
        [125]="F115W"
        [164]="F150W"
        [214]="F210M"
)

declare -a longfilters=(
        [378]="F360M"
        [476]="F480M"
)

# 125 164 214
for lam in "${!shortfilters[@]}"
do
    for dir in 10jup50au 5jup50au 1jup50au 1sat50au 
    do
        for casedir in 5 0 1 2
        do
            for inc in 0 30 60
            do
                python fits_to_adu.py ${dir}/${casedir}/npix84lam${lam}inc${inc}.fits coron
                python nircam_sim.py ${dir}/${casedir}/npix84lam${lam}inc${inc}ADU_coron.fits pipl
                mv image${lam}_sim/jw00005001001_01101_00001_nrcb2_linear_cal.fits image${lam}_sim/${dir}/${casedir}/${shortfilters[$lam]}_lam${lam}inc${inc}_coron_cal.fits
            done
        done
    done

    for dir in 10jup30au 5jup30au 1jup30au
    do
        for casedir in 5 0 1 2
        do
            for inc in 0 30 60
            do
                python fits_to_adu.py ${dir}/${casedir}/npix48lam${lam}inc${inc}.fits coron
                python nircam_sim.py ${dir}/${casedir}/npix48lam${lam}inc${inc}ADU_coron.fits pipl
                mv image${lam}_sim/jw00005001001_01101_00001_nrcb2_linear_cal.fits image${lam}_sim/${dir}/${casedir}/${shortfilters[$lam]}_lam${lam}inc${inc}_coron_cal.fits
            done
        done
    done
done


# 378 476
for lam in "${!longfilters[@]}"
do
    for dir in 10jup50au 5jup50au 1jup50au 1sat50au 
    do
        for casedir in 5 0 1 2
        do
            for inc in 0 30 60
            do
                python fits_to_adu.py ${dir}/${casedir}/npix41lam${lam}inc${inc}.fits coron
                python nircam_sim.py ${dir}/${casedir}/npix41lam${lam}inc${inc}ADU_coron.fits pipl
                mv image${lam}_sim/jw00005001001_01101_00001_nrcb5_linear_cal.fits image${lam}_sim/${dir}/${casedir}/${longfilters[$lam]}_lam${lam}inc${inc}_coron_cal.fits
            done
        done
    done

    for dir in 10jup30au 5jup30au 1jup30au
    do
        for casedir in 5 0 1 2
        do
            for inc in 0 30 60
            do
                python fits_to_adu.py ${dir}/${casedir}/npix23lam${lam}inc${inc}.fits coron
                python nircam_sim.py ${dir}/${casedir}/npix23lam${lam}inc${inc}ADU_coron.fits pipl
                mv image${lam}_sim/jw00005001001_01101_00001_nrcb5_linear_cal.fits image${lam}_sim/${dir}/${casedir}/${longfilters[$lam]}_lam${lam}inc${inc}_coron_cal.fits
            done
        done
    done
done
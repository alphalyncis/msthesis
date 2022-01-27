#!/bin/bash

declare -a filters=(
        [277]="F277W"
        [378]="F380M"
        [430]="F430M"
        [476]="F480M"
)

# # ref stars run
# for lam in "${!filters[@]}"
# do
#     for dir in 10jup50au
#     do
#         for casedir in 5
#         do
#             for inc in 0
#             do
#                 python fits_to_adu.py ${dir}/${casedir}/npix41lam${lam}inc${inc}.fits staronly
#                 python nirissami_sim.py ${dir}/${casedir}/npix41lam${lam}inc${inc}ADU_staronly.fits piplref
#                 mv image${lam}_sim/jw00005001001_01101_00005_nis_linear_cal.fits image${lam}_sim/refstar_nis_linear_cal.fits
#             done
#         done
#     done

#     for dir in 10jup30au
#     do
#         for casedir in 5
#         do
#             for inc in 0
#             do
#                 python fits_to_adu.py ${dir}/${casedir}/npix23lam${lam}inc${inc}.fits staronly
#                 python nirissami_sim.py ${dir}/${casedir}/npix23lam${lam}inc${inc}ADU_staronly.fits piplref
#                 mv image${lam}_sim/jw00005001001_01101_00005_nis_linear_cal.fits image${lam}_sim/refstar_nis_linear_cal.fits
#             done
#         done
#     done
# done


# science run
for lam in "${!filters[@]}"
do
    for dir in 10jup50au 5jup50au 1jup50au 1sat50au
    do
        for casedir in 5 0 1 2 
        do
            for inc in 0 30 60
            do
                python fits_to_adu.py ${dir}/${casedir}/npix41lam${lam}inc${inc}.fits
                python nirissami_sim.py ${dir}/${casedir}/npix41lam${lam}inc${inc}ADU.fits pipl
                mv image${lam}_sim/jw00005001001_01101_00005_nis_linear_cal.fits image${lam}_sim/${dir}/${casedir}/${filters[$lam]}_lam${lam}inc${inc}_ami_cal.fits
                mv image${lam}_sim/combined_aminorm.fits image${lam}_sim/${dir}/${casedir}/${filters[$lam]}_lam${lam}inc${inc}_amiref_cal.fits
            done
        done
    done

    for dir in 10jup30au 5jup30au 1jup30au
    do
        for casedir in 5 0 1 2
        do
            for inc in 0 30 60
            do
                python fits_to_adu.py ${dir}/${casedir}/npix23lam${lam}inc${inc}.fits
                python nirissami_sim.py ${dir}/${casedir}/npix23lam${lam}inc${inc}ADU.fits pipl
                mv image${lam}_sim/jw00005001001_01101_00005_nis_linear_cal.fits image${lam}_sim/${dir}/${casedir}/${filters[$lam]}_lam${lam}inc${inc}_ami_cal.fits
                mv image${lam}_sim/combined_aminorm.fits image${lam}_sim/${dir}/${casedir}/${filters[$lam]}_lam${lam}inc${inc}_amiref_cal.fits
            done
        done
    done
done
#!/bin/bash

# nrc/nrs/miri: use photometry_pipl.py + ${casedir}
# micado/metis: use photometry.py
# micado/metis cutouts: use photometry_cutouts.py + ${casedir}


# nrc miri (pipl)
for instru in miri
do
    for inc in 0 30 60
    do
        for casedir in 5 0 1 2
        do
            python photometry_jwst.py ${instru} ${casedir} ${inc}
        done
    done
done

# nrs ami
for instru in nrsami
do
    for inc in 0 30 60
    do
        for casedir in 5 0 1 2
        do
            python photometry_jwst.py ${instru} ${casedir} ${inc}
        done
    done
done


# metis micado
for instru in metis
do
    for inc in 0 30 60
    do
        for casedir in 5 0 1 2
        do
            python photometry_elt.py ${instru} ${casedir} ${inc}
        done
    done
done

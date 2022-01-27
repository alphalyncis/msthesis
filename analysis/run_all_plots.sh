# #!/bin/bash

# nrc
for casedir in 5
do
    for sep in 30 50
    do
        for inc in 0 30 60
        do
            python plot_cutouts.py nrc ${casedir} ${inc} ${sep} 1 5 nosub
        done
    done
done

# nrs
for casedir in 5
do
    for sep in 30 50
    do
        for inc in 0 30 60
        do
            python plot_cutouts.py nrsami ${casedir} ${inc} ${sep} 0.001 0.2 nosub
        done
    done
done

# miri
for casedir in 5
do
    for sep in 30 50
    do
        for inc in 0 30 60
        do
            python plot_cutouts.py miri ${casedir} ${inc} ${sep} 1 4 nosub
        done
    done
done

# micado
for casedir in 5
do
    for sep in 30 50
    do
        for inc in 0 30 60
        do
            python plot_cutouts.py micado ${casedir} ${inc} ${sep} 4 8
        done
    done
done

# # metis
for casedir in 5
do
    for sep in 30 50
    do
        for inc in 0 30 60
        do
            python plot_cutouts.py metis ${casedir} ${inc} ${sep} 5 8
        done
    done
done


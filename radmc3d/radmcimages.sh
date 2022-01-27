#!/bin/bash

seeds=('s/-27437/-333/g' 's/-333/-93846/g' 's/-93846/-11111/g' 's/-11111/-8112/g' 's/-8112/-27437/g')
nthread=15
declare -a lams50=(
        [125]="1733 84 41"
        [164]="1733 84 41"
        [214]="1733 84 41"
        [277]="41"
        [378]="520 41"
        [430]="41"
        [476]="520 41"
        [1020]="371"
)

declare -a lams30=(
        [125]="1000 48 23"
        [164]="1000 48 23"
        [214]="1000 48 23"
        [277]="23"
        [378]="300 23"
        [430]="23"
        [476]="300 23"
        [1020]="214"
)

declare -a lamscube=(
        [1020]="9. 11."
        [1500]="13. 17.8"
        [2100]="17.5 25"
)



##### non-cube images #####

for lam in "${!lams50[@]}"
do 
    ## 1. repeat 50au ###
    for dir in 10jup50au 5jup50au
    do
        cd ${dir}
        for casedir in 5 0 1 2
        do
            cd ${casedir}
            read -a npixlist <<<"${lams50[$lam]}"
            for n in "${npixlist[@]}"
            do
                for i in 1 2 3 4 5
                do
                    echo ${seeds[i-1]}
                    for a in 0 30 60
                    do
                        wl=`echo $lam/100|bc -l`
                        echo "radmc3d image npix ${n} sizeau 250 incl ${a} lambda ${wl} fluxcons setthreads ${nthread}"
                        echo "$casedir/npix${n}lam${lam}inc${a}_0${i}.fits"
                        cp image.out npix${n}lam${lam}inc${a}_0${i}.out
                        python tofits.py npix${n}lam${lam}inc${a}_0${i}.out
                        rm npix${n}lam${lam}inc${a}_0${i}.out
                    done
                    sed ${seeds[i-1]} radmc3d.inp > idf
                    mv idf radmc3d.inp
                done
            done
            cd ..
        done
        cd ..
    done

    ### 2. non-repeat 50au ###
    for dir in 1jup50au 1sat50au 10jup50au 5jup50au
    do
        cd ${dir}
        for casedir in 0 5 1 2
        do
            cd ${casedir}
            read -a npixlist <<<"${lams50[$lam]}"
            for n in "${npixlist[@]}" 
            do
                for a in 0 30 60
                do
                    wl=`echo $lam/100|bc -l`
                    radmc3d image npix ${n} sizeau 250 incl ${a} lambda ${wl} fluxcons setthreads ${nthread}
                    cp image.out npix${n}lam${lam}inc${a}.out
                    python tofits.py npix${n}lam${lam}inc${a}.out
                    rm npix${n}lam${lam}inc${a}.out
                done
            done
            cd ..
        done
        cd ..    
    done
done


for lam in "${!lams30[@]}"
do
    ### 3. repeat 30au ###
    for dir in 10jup30au 5jup30au
    do
        cd ${dir}
        for casedir in 0 5 1 2
        do
            cd ${casedir}
            read -a npixlist <<<"${lams30[$lam]}"
            for n in "${npixlist[@]}"
            do
                for i in 1 2 3 4 5
                do
                    echo ${seeds[i-1]}
                    for a in 0 30 60
                    do
                        wl=`echo $lam/100|bc -l`
                        radmc3d image npix ${n} sizeau 150 incl ${a} lambda ${wl} fluxcons setthreads ${nthread}
                        cp image.out npix${n}lam${lam}inc${a}_0${i}.out
                        python tofits.py npix${n}lam${lam}inc${a}_0${i}.out
                        rm npix${n}lam${lam}inc${a}_0${i}.out
                    done
                    sed ${seeds[i-1]} radmc3d.inp > idf
                    mv idf radmc3d.inp
                done
            done
            cd ..
        done
        cd ..
    done

    ### 4. non-repeat 30au ###
    for dir in 1jup30au
    do
        cd ${dir}
        for casedir in 0 5 1 2
        do
            cd ${casedir}
            read -a npixlist <<<"${lams30[$lam]}"
            for n in "${npixlist[@]}"
            do
                for a in 0 30 60
                do
                    wl=`echo $lam/100|bc -l`
                    radmc3d image npix ${n} sizeau 150 incl ${a} lambda ${wl} fluxcons setthreads ${nthread}
                    cp image.out npix${n}lam${lam}inc${a}.out
                    python tofits.py npix${n}lam${lam}inc${a}.out
                    rm npix${n}lam${lam}inc${a}.out
                done
            done
            cd ..
        done
        cd ..    
    done
done

# remember to run ""python avgfits.py ${lam}"" afterwards





##### cube images #####


for lam in "${!lamscube[@]}"
do
    wl="${lamscube[$lam]}"

    ### 1. repeat 50au ###
    for dir in 10jup50au 5jup50au
    do
        cd ${dir}
        for casedir in 5 0 1 2
        do
            cd ${casedir}
            for n in 23
            do
                for i in 1 2 3 4 5
                do
                    echo ${seeds[i-1]}
                    for a in 0 30 60
                    do
                        radmc3d image npix ${n} sizeau 250 incl ${a} lambdarange ${wl} nlam 5 fluxcons setthreads ${nthread}
                        cp image.out cube5npix${n}lam${lam}inc${a}_0${i}.out
                        python tofits_mirisim.py cube5npix${n}lam${lam}inc${a}_0${i}.out
                        rm cube5npix${n}lam${lam}inc${a}_0${i}.out
                    done
                    sed ${seeds[i-1]} radmc3d.inp > idf
                    mv idf radmc3d.inp
                done
            done
            cd ..
        done
        cd ..
    done

    ### 2. non-repeat 50au ###
    for dir in 1jup50au 1sat50au
    do
        cd ${dir}
        for casedir in 5 0 1 2
        do
            cd ${casedir}
            for n in 23
            do
                for a in 0 30 60
                do
                    radmc3d image npix ${n} sizeau 250 incl ${a} lambdarange ${wl} nlam 5 fluxcons setthreads ${nthread}
                    cp image.out cube5npix${n}lam${lam}inc${a}.out
                    python tofits_mirisim.py cube5npix${n}lam${lam}inc${a}.out
                    rm cube5npix${n}lam${lam}inc${a}.out
                done
            done
            cd ..
        done
        cd ..
    done

    ### 3. repeat 30au ###
    for dir in 10jup30au 5jup30au
    do
        cd ${dir}
        for casedir in 5 0 1 2
        do
            cd ${casedir}
            for n in 13
            do
                for i in 1 2 3 4 5
                do
                    echo ${seeds[i-1]}
                    for a in 0 30 60
                    do
                        radmc3d image npix ${n} sizeau 150 incl ${a} lambdarange ${wl} nlam 5 fluxcons setthreads ${nthread}
                        cp image.out cube5npix${n}lam${lam}inc${a}_0${i}.out
                        python tofits_mirisim.py cube5npix${n}lam${lam}inc${a}_0${i}.out
                        rm cube5npix${n}lam${lam}inc${a}_0${i}.out
                    done
                    sed ${seeds[i-1]} radmc3d.inp > idf
                    mv idf radmc3d.inp
                done
            done
            cd ..
        done
        cd ..
    done

    ### 4. non-repeat 30au ###
    for dir in 1jup30au
    do
        cd ${dir}
        for casedir in 5 0 1 2
        do
            cd ${casedir}
            for n in 13
            do
                for a in 0 30 60
                do
                    radmc3d image npix ${n} sizeau 150 incl ${a} lambdarange ${wl} nlam 5 fluxcons setthreads ${nthread}
                    cp image.out cube5npix${n}lam${lam}inc${a}.out
                    python tofits_mirisim.py cube5npix${n}lam${lam}inc${a}.out
                    rm cube5npix${n}lam${lam}inc${a}.out
                done
            done
            cd ..
        done
        cd ..
    done
done

# python avgfits.py ${lam} cube
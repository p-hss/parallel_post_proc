reset

set terminal epslatex color size 6,4.8 

set xrange[:]
set yrange[:]


set loadpath '../../../.config'
load 'spectral.pal'

set ylabel 'Histogram'
#set key left top 
set key below
#set key outside center right

mhd=1 
sim="R01 R02 R03 R04 R05 R06 R07 R08"

time="20"
files_1 = "../data/sim_"
files_2 = "/angle_histo_"
files_3 = "_tkolmo_".time.".dat"

if(mhd == 0){
    set output "figures/histograms/hd_angle_histo_t".time.".tex"
    set xlabel '$|\cos(\Gamma(\va{T}_i,\va{l}))|$'
    plot for [i=1:words(sim)]\
        files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:2 w l \
        title word(sim,i).': $H(|\cos(\measuredangle(\va{T}_1,\va{l}))|)$'  ls i+9,\
        for [i=1:words(sim)]\
        files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:3 w l \
        title word(sim,i).': $H(|\cos(\measuredangle(\va{T}_2,\va{l}))|)$'  ls i+19,\
        for [i=1:words(sim)]\
        files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:4 w l \
        title word(sim,i).': $H(|\cos(\measuredangle(\va{T}_3,\va{l}))|)$'  ls i+29,\
        for [i=1:words(sim)]\
        files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:5 w l \
        title word(sim,i).': $H(|\cos(\measuredangle(\va{\omega},\va{l}))|)$'  ls i+39,\
}

if(mhd == 1){
    set output "figures/histograms/mhd_angle_histo_t".time.".tex"
    set xlabel '$|\cos(\Gamma(\va{T}_i,\va{l}))|$'
    plot for [i=1:words(sim)]\
        files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:2 w l \
        title word(sim,i).': $H(|\cos(\measuredangle(\va{T}_1,\va{l}))|)$'  ls i+9,\
        for [i=1:words(sim)]\
        files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:6 w l \
        title word(sim,i).': $H(|\cos(\measuredangle(\delta \va{b},\va{l}))|)$'  ls i+19,\
}

        #for [i=1:words(sim)]\
        #files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:3 w l \
        #title word(sim,i).': $H(|\cos(\measuredangle(\va{T}_2,\va{l}))|)$'  ls i+19,\
        #for [i=1:words(sim)]\
        #files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:4 w l \
        #title word(sim,i).': $H(|\cos(\measuredangle(\va{T}_3,\va{l}))|)$'  ls i+29,\
        #for [i=1:words(sim)]\
        #files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:5 w l \
        #title word(sim,i).': $H(|\cos(\measuredangle(\va{\omega},\va{l}))|)$'  ls i+39,\
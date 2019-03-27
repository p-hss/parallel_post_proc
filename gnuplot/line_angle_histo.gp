reset

set terminal epslatex color size 6,4 

set xrange[:]
set yrange[:]


set loadpath '../../gnuplot_palettes/'
load 'spectral_poster.pal'

set ylabel 'Probability density function'
#set key left top 
set key top left Left reverse
#set key outside center right

mhd=1 
#sim="R01 R02 R03 R04 R05 R06"
sim="R01"

time="20"
f_1 = "../data/sim_"
f_2 = "/angle_histo_"
f_3 = "_tkolmo_".time.".dat"

g_1 = "../data/sim_"
g_2 = "/mhd_angle_histo_"
g_3 = "_tkolmo_".time.".dat"


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
    #set xlabel '$|\va{\hat{T}}_i \cdot \va{\hat{l}}|$'
    set xlabel '$|\va{\hat{l}} \cdot \va{\hat{x}}| \qquad \va{\hat{x}} = \va{\hat{T}}_i, \, \va{\hat{\omega}}, \, \va{\hat{b}} $'

    plot for [i=1:words(sim)]\
        f_1.word(sim,i).f_2.word(sim,i).f_3 u 1:2 w l \
        title '$|\va{\hat{l}} \cdot \va{\hat{T}}_1|$'  ls 100,\
        f_1.word(sim,i).f_2.word(sim,i).f_3 u 1:3 w l \
        title '$|\va{\hat{l}} \cdot \va{\hat{T}}_2|$'  ls 200,\
        f_1.word(sim,i).f_2.word(sim,i).f_3 u 1:5 w l \
        title '$|\va{\hat{l}} \cdot \va{\hat{\omega}}|$'  ls 300,\
        g_1.word(sim,i).g_2.word(sim,i).g_3 u 1:2 w l \
        title '$|\va{\hat{l}} \cdot \va{\hat{b}}|$'  ls 400,\
}

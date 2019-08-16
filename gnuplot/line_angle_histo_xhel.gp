reset

set terminal epslatex color size 6,4 

set xrange[:]
set yrange[:]

set loadpath '../../gnuplot_palettes/'
#load 'spectral_poster.pal'
load 'contrast.pal'

set ylabel 'Probability density function'
#set key left top 
#set key top left Left reverse
set key below Left reverse
set key spacing 1.15
#set key outside center right

sim="Z41 Z43 Z46"
names="$\\sigma_C=1$ $\\sigma_C=0.5$ $\\sigma_C=0$"

time="20"

f_1 = "../data/sim_"
f_2 = "/angle_histo_"
f_3 = "_tkolmo_".time.".dat"

g_1 = "../data/sim_"
g_2 = "/mhd_angle_histo_"
g_3 = "_tkolmo_".time.".dat"

    #plot file1 u 1:3 w l \
    #     title sim1.': $H(|\cos(\measuredangle(\delta \va{b},\va{T_1}))|)$' ls 18, \
    #     file1 u 1:4 w l \

set output "figures/histograms/xhel_angle_histo_t".time."_256.tex"
set xlabel 'angle'
plot for [i=1:words(sim)]\
    f_1.word(sim,i).f_2.word(sim,i).f_3 u 1:2 w l \
    title word(names,i).': $H(|\va{\hat{l}} \cdot \va{\hat{T}}_1|)$'  ls i+9,\
    for [i=1:words(sim)]\
    g_1.word(sim,i).g_2.word(sim,i).g_3 u 1:3 w l \
    title word(names,i).': H$(|\va{\hat{b}} \cdot \va{\hat{T}}_1|)$'  ls i+29,\
#plot for [i=1:words(sim)]\
#    f_1.word(sim,i).f_2.word(sim,i).f_3 u 1:2 w l \
#    title word(names,i).': $H(|\va{\hat{l}} \cdot \va{\hat{T}}_1|)$'  ls i+9,\
#    for [i=1:words(sim)]\
#    g_1.word(sim,i).g_2.word(sim,i).g_3 u 1:2 w l \
#    title word(names,i).': H$(|\va{\hat{l}} \cdot \va{\hat{b}}|)$'  ls i+19,\
#    for [i=1:words(sim)]\
#    g_1.word(sim,i).g_2.word(sim,i).g_3 u 1:3 w l \
#    title word(names,i).': H$(|\va{\hat{T}}_1 \cdot \va{\hat{b}}|)$'  ls i+29,\

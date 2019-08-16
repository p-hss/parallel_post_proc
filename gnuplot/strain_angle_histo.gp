reset

set terminal epslatex color size 6,4 

set xrange[:]
set yrange[:]

set loadpath '../../gnuplot_palettes/'
load 'spectral_poster.pal'
#load 'contrast.pal'

set ylabel 'Probability density function'
set xlabel 'Angle'
set key below

sim="M04 H00"
names="MHD HD"
time="20"

s1="100 1000"
s2="200 2000"
s3="300 3000"
s4="400 4000"
s5="500 5000"

f_1 = "../data/sim_"
f_2 = "/angle_histo_"
f_3 = "_tkolmo_".time.".dat"

g_1 = "../data/sim_"
g_2 = "/mhd_angle_histo_"
g_3 = "_tkolmo_".time.".dat"

set output "figures/histograms/strain_angle_histo_t20_256.tex"
plot for [i=1:words(sim)]\
    f_1.word(sim,i).f_2.word(sim,i).f_3 u 1:6 w l \
    title word(names,i).': $|\va{\hat{T}}_1 \cdot \va{\hat{\omega}}|$' ls word(s1,i),\
    for [i=1:words(sim)]\
    f_1.word(sim,i).f_2.word(sim,i).f_3 u 1:7 w l \
    title word(names,i).': $|\va{\hat{T}}_2 \cdot \va{\hat{\omega}}|$' ls word(s2,i),\
    g_1.word(sim,1).g_2.word(sim,1).g_3 u 1:3 w l \
    title word(names,1).': $|\va{\hat{T}}_1 \cdot \va{b}|$' ls word(s3,1),\
    g_1.word(sim,1).g_2.word(sim,1).g_3 u 1:4 w l \
    title word(names,1).': $|\va{\hat{T}}_2 \cdot \va{b}|$' ls word(s4,1)

    #plot file1 u 1:3 w l \
    #     title sim1.': $H(|\cos(\measuredangle(\delta \va{b},\va{T_1}))|)$' ls 18, \
    #     file1 u 1:4 w l \
    #     title sim1.': $H(|\cos(\measuredangle(\delta \va{b},\va{T_2}))|)$' ls 28, \
    #     file1 u 1:5 w l \
    #     title sim1.': $H(|\cos(\measuredangle(\delta \va{b},\va{T_3}))|)$' ls 38, \
    #     file1 u 1:6 w l \
    #     title sim1.': $H(|\cos(\measuredangle(\delta \va{b},\va{\omega}))|)$' ls 48, \
    #     file2 u 1:3 w l \
    #     title sim2.': $H(|\cos(\measuredangle(\delta \va{b},\va{T_1}))|)$' ls 17, \
    #     file2 u 1:4 w l \
    #     title sim2.': $H(|\cos(\measuredangle(\delta \va{b},\va{T_2}))|)$' ls 27, \
    #     file2 u 1:5 w l \
    #     title sim2.': $H(|\cos(\measuredangle(\delta \va{b},\va{T_3}))|)$' ls 37, \
    #     file2 u 1:6 w l \
    #     title sim2.': $H(|\cos(\measuredangle(\delta \va{b},\va{\omega}))|)$' ls 47, \


    #set xlabel '$|\cos(\Gamma(\va{T}_i,\va{\omega}))|$'
    #set output "figures/histograms/mhd_strain_vorticity_angle_histo_t20.tex"

    #plot file1 u 1:7 w l\
    #     title sim1.': $H(|\cos(\measuredangle(\va{T_1},\va{\omega}))|)$' ls 18, \
    #     file1 u 1:8 w l \
    #     title sim1.': $H(|\cos(\measuredangle(\va{T_2},\va{\omega}))|)$' ls 28, \
    #     file1 u 1:9 w l \
    #     title sim1.': $H(|\cos(\measuredangle(\va{T_3},\va{\omega}))|)$' ls 38, \
    #     file2 u 1:7 w l\
    #     title sim2.': $H(|\cos(\measuredangle(\va{T_1},\va{\omega}))|)$' ls 17, \
    #     file2 u 1:8 w l \
    #     title sim2.': $H(|\cos(\measuredangle(\va{T_2},\va{\omega}))|)$' ls 27, \
    #     file2 u 1:9 w l \
    #     title sim2.': $H(|\cos(\measuredangle(\va{T_3},\va{\omega}))|)$' ls 37, \
reset

set terminal epslatex color size 6,4 

set xrange[:]
set yrange[:]

set loadpath '../../gnuplot_palettes/'
load 'spectral_poster.pal'
#load 'contrast.pal'

set ylabel 'Probability density function'
#set key left top 
#set key top left Left reverse
set key below
#set key spacing 1.5
#set key outside center right

mhd=1 
#sim="R01 R02 R03 R04 R05 R06"
sim="X40 X46"
names="XHEL1 XHEL6"

s1="100 1000"

time="20"
f_1 = "../data/sim_"
f_2 = "/angle_histo_"
f_3 = "_tkolmo_".time.".dat"

g_1 = "../data/sim_"
g_2 = "/mhd_angle_histo_"
g_3 = "_tkolmo_".time.".dat"

s1="100 1000"
s2="200 2000"
s3="300 3000"
s4="400 4000"
s5="500 5000"

set output "figures/histograms/xhel_strain_angle_histo_t".time."_256.tex"
set xlabel 'Alignment'

plot for [i=1:words(sim)]\
    f_1.word(sim,i).f_2.word(sim,i).f_3 u 1:2 w l \
    title word(names,i).': $|\va{\hat{l}} \cdot \va{\hat{T}}_1|$' ls word(s1,i),\
    for [i=1:words(sim)]\
    f_1.word(sim,i).f_2.word(sim,i).f_3 u 1:3 w l \
    title word(names,i).': $|\va{\hat{l}} \cdot \va{\hat{T}}_2|$' ls word(s2,i),\
    for [i=1:words(sim)]\
    f_1.word(sim,i).f_2.word(sim,i).f_3 u 1:4 w l \
    title word(names,i).': $|\va{\hat{l}} \cdot \va{\hat{T}}_3|$' ls word(s4,i),\

set output "figures/histograms/xhel_fields_angle_histo_t".time."_256.tex"
plot for [i=1:words(sim)]\
    f_1.word(sim,i).f_2.word(sim,i).f_3 u 1:5 w l \
    title word(names,i).': $|\va{\hat{l}} \cdot \va{\hat{\omega}}|$' ls word(s1,i),\
    for [i=1:words(sim)]\
    g_1.word(sim,i).g_2.word(sim,i).g_3 u 1:2 w l \
    title word(names,i).': $|\va{\hat{l}} \cdot \va{\hat{b}}|$'  ls word(s2,i),\
    for [i=1:words(sim)]\
    g_1.word(sim,i).g_2.word(sim,i).g_3 u 1:12 w l \
    title word(names,i).': $|\va{\hat{l}} \cdot \va{\hat{v}}|$'  ls word(s3,i)

set output "figures/histograms/xhel_b_strain_angle_histo_t".time."_256.tex"
plot for [i=1:words(sim)]\
    g_1.word(sim,i).g_2.word(sim,i).g_3 u 1:3 w l \
    title word(names,i).': $|\va{\hat{b}} \cdot \va{\hat{T}}_1|$'  ls word(s1,i),\
    for [i=1:words(sim)]\
    g_1.word(sim,i).g_2.word(sim,i).g_3 u 1:4 w l \
    title word(names,i).': $|\va{\hat{b}} \cdot \va{\hat{T}}_2|$'  ls word(s2,i),\
    for [i=1:words(sim)]\
    g_1.word(sim,i).g_2.word(sim,i).g_3 u 1:5 w l \
    title word(names,i).': $|\va{\hat{b}} \cdot \va{\hat{T}}_3|$'  ls word(s3,i)

set logscale y
set yrange[0.1:160]
set output "figures/histograms/xhel_mixed_fields_angle_histo_t".time."_256.tex"
plot for [i=1:words(sim)]\
    g_1.word(sim,i).g_2.word(sim,i).g_3 u 1:10 w l \
    title word(names,i).': $|\va{\hat{v}} \cdot \va{\hat{b}}|$'  ls word(s1,i),\
    for [i=1:words(sim)]\
    g_1.word(sim,i).g_2.word(sim,i).g_3 u 1:11 w l \
    title word(names,i).': $|\va{\hat{v}} \cdot \va{\hat{T}}_1|$'  ls word(s3,i)
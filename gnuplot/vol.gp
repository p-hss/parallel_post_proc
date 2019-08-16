reset

set terminal epslatex color size 6,4 

set xrange[0:30]
set yrange[:0]

set loadpath '../../gnuplot_palettes/'
load 'spectral_poster.pal'
#load 'contrast.pal'

set xlabel '$t/ \tau_{\eta}$'
set ylabel '$\langle \Theta \rangle, \quad \langle  \Phi \rangle$'

set key right bottom 
unset grid

s1="100 1000"
s2="200 2000"

sim="M04 H00"
names="MHD HD"

time = "20"
files_1 = "../data/sim_"
files_2 = "/line_evo_"
files_3 = "/surf_evo_"
files_4 = ".dat"

set output "figures/vol_evo_".time."_256.tex"
plot for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_4 u 1:4 w l \
     title word(names,i).': $\langle \Theta \rangle$ ' ls word(s1,i),\
     for [i=1:words(sim)]\
     files_1.word(sim,i).files_3.word(sim,i).files_4 u 1:4 w l \
     title word(names,i).': $\langle \Phi \rangle$ ' ls word(s2,i),\


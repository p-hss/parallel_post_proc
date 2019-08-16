reset

set terminal epslatex color size 6,4 

set loadpath '../../gnuplot_palettes/'
#load 'contrast.pal'
load 'spectral_poster.pal'

set yrange[0:0.2]
set xrange[0:30]

set xlabel '$t/ \tau_{\eta}$'
set ylabel '$\langle \xi \rangle, \quad \textrm{var}(\xi)$ '

set key
set key spacing 1.25

sim = "M04 H00"
names = "MHD HD"

files_1 = "../data/sim_"
files_2 = "/surf_evo_"
files_3 = ".dat"

s1="100 1000"
s2="200 2000"

set output "figures/surface_evo_256.tex"
plot for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:2 w l \
     title word(names,i).' $\langle \xi \rangle$' ls word(s1,i),\
     for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:3 w l \
     title word(names,i).' var$(\xi)$' ls word(s2,i),\


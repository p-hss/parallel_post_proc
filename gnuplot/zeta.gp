reset

set terminal epslatex dashed color size 6,4 

set xrange[0:30]
set yrange[0:0.18]
#set y2range[0.7:1]
#set ytics 0.05
#set y2tics 0.1

set xlabel '$t/ \tau_{\eta}$'
set ylabel '$\langle \zeta \rangle, \quad \textrm{var}(\zeta) $'
#set y2label '$\langle \sphericalangle(\va{l} ,\va{T}_1) \rangle $ '

set loadpath '../../gnuplot_palettes/'
#load 'contrast.pal'
load 'spectral_poster.pal'

set key
set key spacing 1.25

sim = "M04 H00"
names = "MHD HD"

s1="100 1000"
s2="200 2000"

files_1 = "../data/sim_"
files_2 = "/line_evo_"
files_3 = ".dat"


set output "figures/line_evo_256.tex"
plot for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:2 w l \
     title word(names,i).' $\langle \zeta \rangle$' ls word(s1,i),\
     for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:3 w l \
     title word(names,i).' var$(\zeta)$' ls word(s2,i),\

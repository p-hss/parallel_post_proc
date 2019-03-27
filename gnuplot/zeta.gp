reset

set terminal epslatex color size 6,4 

set xrange[0:30]
set yrange[0:0.20]
set y2range[0.7:1]
set ytics 0.05
set y2tics 0.1

set xlabel '$t/ \tau_{\eta}$'
set ylabel '$\langle \zeta \rangle$ '
set y2label '$\langle \sphericalangle(\va{l} ,\va{T}_1) \rangle $ '
set ytics nomirror

set loadpath '../../gnuplot_palettes/'
load 'spectral_poster.pal'

set key
sim = "Z01"

files_1 = "../data/sim_"
files_2 = "/line_evo_"
files_3 = ".dat"

f_1 = "../data/sim_"
f_2 = "/gamma_"
f_3 = ".dat"

set output "figures/mhd_line_evo_256.tex"
plot for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:2 w l axis x1y1 \
     title '$\langle \zeta \rangle$' ls 100,\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:3 w l axis x1y1 \
     title 'var$(\zeta)$' ls 200,\
     for [i=1:words(sim)]\
          f_1.word(sim,i).f_2.word(sim,i).f_3 u 1:($2) w l \
          title '$\langle \sphericalangle(\va{l},\va{T}_1) \rangle$' ls 300 axis x1y2

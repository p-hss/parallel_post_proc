reset

set terminal epslatex color size 6,4 

set xrange[0:30]
set yrange[:0]

set loadpath '../../gnuplot_palettes/'
load 'spectral_poster.pal'

set xlabel '$t/ \tau_{\eta}$'
set ylabel '$\langle \Theta \rangle \quad \
\langle  \Phi \rangle$'

set key right bottom 
unset grid

sim = "M"

time = "20"
files_1 = "../data/sim_"
files_2 = "/line_evo_"
files_3 = "/surf_evo_"
files_4 = ".dat"

set output "figures/histograms/vol_evo".time."_256.pdf"
plot for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_4 u 1:4 w l \
     title '$\langle \Theta \rangle$ ' ls i+9,\
     for [i=1:words(sim)]\
     files_1.word(sim,i).files_3.word(sim,i).files_4 u 1:4 w l \
     title '$\langle \Phi \rangle$ ' ls i+19,\


#if(mhd == 0){
#
#    load 'blues.pal'
#    set output "figures/hd_vol_evo.tex"
#    dir = "../data/sim_"
#    file = ".dat"
#    name1 = "/line_evo_"
#    name2 = "/surf_evo_"
#
#    sim1 = "C05"
#    file11 = dir.sim1.name1.sim1.file
#    file12 = dir.sim1.name2.sim1.file
#    sim2 = "C07"
#    file21 = dir.sim2.name1.sim2.file
#    file22 = dir.sim2.name2.sim2.file
#
#    plot file11 u 1:4 w l  title sim1.': $\langle \Theta \rangle$' ls 18 ,\
#         file12 u 1:4 w l  title sim1.': $\langle  \Phi \rangle$' ls 38,\
#         file21 u 1:4 w l  title sim2.': $\langle \Theta \rangle$' ls 17 ,\
#         file22 u 1:4 w l  title sim2.': $\langle  \Phi \rangle$' ls 37
#}
#
#if(mhd == 1){
#
#    load 'reds.pal'
#    set output "figures/mhd_vol_evo.tex"
#    dir = "../data/sim_"
#    file = ".dat"
#    name1 = "/line_evo_"
#    name2 = "/surf_evo_"
#
#    sim1 = "Q02"
#    file11 = dir.sim1.name1.sim1.file
#    file12 = dir.sim1.name2.sim1.file
#    sim2 = "Q10"
#    file21 = dir.sim2.name1.sim2.file
#    file22 = dir.sim2.name2.sim2.file
#
#    plot file11 u 1:4 w l  title sim1.': $\langle \Theta \rangle$' ls 18,\
#         file12 u 1:4 w l  title sim1.': $\langle  \Phi \rangle$' ls 38,\
#         file21 u 1:4 w l  title sim2.': $\langle \Theta \rangle$' ls 17,\
#         file22 u 1:4 w l  title sim2.': $\langle  \Phi \rangle$' ls 37 
#}


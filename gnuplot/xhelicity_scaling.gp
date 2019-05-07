reset

set terminal epslatex color size 6,4 

set xrange[:]
set yrange[:]

set loadpath '../../gnuplot_palettes/'
load 'spectral_poster.pal'

set xlabel '$\sigma_c$'
#set ylabel '$\zeta/\zeta_{\textrm{min}}$, $\sphericalangle(\va{l}, \va{T_1})/\sphericalangle(\va{l}, \va{T_1})_{\textrm{max}}$, $\sphericalangle(\va{l}, \va{T_1})/\sphericalangle(\va{l}, \va{b})_{\textrm{max}}$'
set ylabel ''
set key left top Left reverse

set arrow from graph 0, first 1 to graph 1, first 1 nohead lw 1.5 lc "gray"
set ytics 0.1
set key spacing 1.5

file = "../data/line_xhel_scaling_256.dat"

set output "figures/line_xhel_Z4_scaling_256.tex"
plot file u 1:2:3 w yerrorbars notitle ls 100,\
    file u 1:2 w l  title '$\overline{\zeta}/\overline{\zeta}_{\sigma_{c}=0}$' ls 100,\
    file u 1:4:5 w yerrorbars notitle  ls 100,\
    file u 1:4 w l title '$\sphericalangle(\va{l}, \va{T_1})/\sphericalangle(\va{l}, \va{T_1})_{\sigma_{c}=0}$'  ls 200,\
    file u 1:6:7 w yerrorbars notitle  ls 100,\
    file u 1:6 w l title '$\sphericalangle(\va{l}, \va{b})/\sphericalangle(\va{l}, \va{b})_{\sigma_{c}=0}$'  ls 300,\

reset

set terminal epslatex color size 6,4 

set xrange[:]
set yrange[:]

set loadpath '../../gnuplot_palettes/'
load 'spectral.pal'

set xlabel '$\sigma_C$'
#set ylabel '$\zeta/\zeta_{\textrm{min}}$, $\sphericalangle(\va{l}, \va{T_1})/\sphericalangle(\va{l}, \va{T_1})_{\textrm{max}}$, $\sphericalangle(\va{l}, \va{T_1})/\sphericalangle(\va{l}, \va{b})_{\textrm{max}}$'
set ylabel ''
set key left top 

set arrow from graph 0, first 1 to graph 1, first 1 nohead lc "black"

file = "../data/line_xhel_scaling.dat"

set output "figures/line_xhel_scaling_128.tex"
plot file u 1:2:3 w yerrorbars notitle ls 100,\
    file u 1:2 w l  title '$\zeta/\zeta_{\sigma_{c,f}=0}$' ls 100,\
    file u 1:4:5 w yerrorbars notitle  ls 100,\
    file u 1:4 w l title '$\sphericalangle(\va{l}, \va{T_1})/\sphericalangle(\va{l}, \va{T_1})_{\textrm{max}}$'  ls 200,\
    file u 1:6:7 w yerrorbars notitle  ls 100,\
    file u 1:6 w l title '$\sphericalangle(\va{l}, \va{T_1})/\sphericalangle(\va{l}, \va{b})_{\textrm{max}}$'  ls 300,\

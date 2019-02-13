reset
set title "Lagrangian particle trajectory"

#set autoscale

#set xrange[:]
#set yrange[1.4:2.6]

unset border
unset xtics
unset ytics
unset ztics
set arrow from 0,0,0 to 0,128,0 nohead
set arrow from 0,0,0 to 0,0,128 nohead
set arrow from 0,0,0 to 128,0,0 nohead

set grid 
#set view 50,30,1,1 
set ticslevel 0

set style line 1 lt 1 lc rgb '#000000' # black 
set style line 2 pt 7 lc rgb '#00008B' ps 0.4  # dark blue
set style line 3 pt 7 lc rgb '#8B0000' ps 1 # dark red 
set style line 4 lt 1 lc rgb '#66CDAA' # grey 
set style line 5 lt 1 lc rgb '#EE82EE' # violet 

set key

set terminal epslatex color size 6,4 
#set output 'figures/lp_pos_trajectory.tex'

filename = "../data/lp_pos_3.dat"

splot filename u 1:2:3 w p title 'particle position' ls 2


reset

#set terminal epslatex color size 6,6 
set terminal png
set yrange[-5:10]
set xrange[0:15]
set loadpath '../../../.config'
load 'spectral.pal'
set palette negative

set size square

set xlabel '$\ln(w_1)$'
set ylabel '$\ln(w_2)$'
unset ztics

unset key
set style line 50 lt 1 lw 1 lc rgb 'black'
set style line 51 lt 3 lw 1 lc rgb 'black'

#p(x)= x
#c(x)= -0.5*x

mhd = 0 

if(mhd == 0){

    sim = "B05"
    dir = "../data/sim_B05"
    set loadpath dir
    name = "/cg_histo_"
    name2 = "/tkolmo_20"
    file = ".dat"
    file1 = name.sim.name2.file

    set output "figures/histograms/hd_cg_histo_t20.png"

    f(x) = x>0.00000001 ? x : 1/0
    set style fill transparent solid 1 noborder
    set style circle radius 0.001
    plot '<sort -g -k3 ../data/sim_B05/cg_histo_'.sim.'_tkolmo_20.dat' u 1:2:(0.10):(f($3))\
    w circles lc palette
    #title '', p(x) w l title "Pancake shape" ls 50, c(x) title "Cigar shape" ls 51

}

if(mhd == 1){

    sim = "Q02"
    dir = "../data/sim_"
    set loadpath dir.sim
    name = "/cg_histo_"
    name2 = "/tkolmo_20"
    file = ".dat"
    file1 = name.sim.name2.file

    set output "figures/histograms/mhd_cg_histo_t20.tex"

    f(x) = x>0.00000000000000001 ? x : 1/0
    set style fill transparent solid 1 noborder
    set style circle radius 0.001
    plot '<sort -g -k3 ../data/sim_Q02/cg_histo_'.sim.'_tkolmo_20.dat' u 1:2:(0.10):(f($3))\
    w circles lc palette
    #title '', p(x) w l title "Pancake shape" ls 50, c(x) title "Cigar shape" ls 51

}

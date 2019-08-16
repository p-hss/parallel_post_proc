reset
set terminal epslatex color size 6,4 
set loadpath '../../gnuplot_palettes/'
load 'spectral_poster.pal'
#load 'contrast.pal'

set key below
set zeroaxis
set xrange[0:16]
set yrange[:]
set xlabel '$t/\tau_{\eta}$'
set ylabel '$\langle \ln(w_1) \rangle, \quad  \langle \ln(w_2) \rangle, \quad \
\langle \ln(w_3) \rangle$'

set bmargin 6

mhd=3 

files_1 = "../data/sim_"
files_2 = "/cg_"
files_3 = ".dat"

sim="M04 H00"
names="MHD HD"

s1="100 1000"
s2="200 2000"
s4="400 4000"

time = "20"
files_1 = "../data/sim_"
files_2 = "/cg_"
files_3 = ".dat"

set output "figures/cg.tex"
plot for [i=1:words(sim)]\
        files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:2 w l \
        title word(names,i).': $\langle \ln(w_1) \rangle$'  ls word(s1,i),\
    for [i=1:words(sim)]\
        files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:3 w l \
        title word(names,i).': $\langle \ln(w_2) \rangle$'  ls word(s2,i),\
    for [i=1:words(sim)]\
        files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:4 w l \
        title word(names,i).': $\langle \ln(w_3) \rangle$'  ls word(s4,i)

if(mhd == 0){
    set output "figures/hd_cg.tex"
    plot for [i=1:words(sim)]\
        files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:2 w l \
        title word(sim,i).': $\langle \ln(w_1) \rangle$'  ls i+9,\
    for [i=1:words(sim)]\
            files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:3 w l \
            title word(sim,i).': $\langle \ln(w_2) \rangle$'  ls i+19,\
    for [i=1:words(sim)]\
            files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:4 w l \
            title word(sim,i).': $\langle \ln(w_3) \rangle$'  ls i+29,\

}
if(mhd == 1){
    set output "figures/mhd_cg.tex"
    plot for [i=1:words(sim)]\
        files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:2 w l \
        title word(sim,i).': $\langle \ln(w_1) \rangle$'  ls i+9,\
    for [i=1:words(sim)]\
            files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:3 w l \
            title word(sim,i).': $\langle \ln(w_2) \rangle$'  ls i+19,\
    for [i=1:words(sim)]\
            files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:4 w l \
            title word(sim,i).': $\langle \ln(w_3) \rangle$'  ls i+29,\
}

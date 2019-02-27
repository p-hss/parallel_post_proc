reset
set terminal epslatex color size 6,4 
set loadpath '../../../.config'
load 'spectral.pal'
set key below
set zeroaxis
set xrange[0:16]
set yrange[:]
set xlabel '$t/\tau_{\eta}$'
set ylabel '$\langle \ln(w_1) \rangle, \quad  \langle \ln(w_2) \rangle, \quad \
\langle \ln(w_3) \rangle$'

#====================================================

mhd=1 
sim="R01 R02 R05 R07 R08"

files_1 = "../data/sim_"
files_2 = "/cg_"
files_3 = ".dat"

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

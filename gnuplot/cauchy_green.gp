reset

set title "Cauchy-Green tensor"

set terminal epslatex color size 6,4 

set xrange[0:16]
set yrange[:]

set loadpath '../../../.config'
load 'blues.pal'

set xlabel '$t/\tau_{\eta}$'
set ylabel '$\langle \ln(w_1) \rangle, \quad  \langle \ln(w_2) \rangle, \quad \
\langle \ln(w_3) \rangle$'

set key left bottom
set zeroaxis

mhd = 0 

if(mhd == 0){

    load 'blues.pal'
    sim = "B05"
    dir = "../data/sim_"
    set loadpath dir.sim
    name = "/cg_"
    file = ".dat"
    file1 = name.sim.file

    set output "figures/hd_cg.tex"

    plot file1 u 1:2 w l title 'B5: $\langle \ln(w_1) \rangle$' ls 18,\
        file1 u 1:3 w l  title 'B5: $\langle \ln(w_2) \rangle$' ls 28,\
        file1 u 1:4 w l  title 'B5: $\langle \ln(w_3) \rangle$' ls 38,\
}

if(mhd == 1){

    load 'reds.pal'
    sim = "m02"
    dir = "../data/sim_"
    set loadpath dir.sim
    name = "/cg_"
    file = ".dat"
    file1 = name.sim.file

    set output "figures/mhd_cg.tex"

    plot file1 u 1:2 w l title 'M2: $\langle \ln(w_1) \rangle$' ls 18,\
        file1 u 1:3 w l  title 'M2: $\langle \ln(w_2) \rangle$' ls 28,\
        file1 u 1:4 w l  title 'M2: $\langle \ln(w_3) \rangle$' ls 38,\
}


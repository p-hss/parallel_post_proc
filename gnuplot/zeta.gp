reset

#set terminal epslatex color size 6,4 
set terminal png

set xrange[0:]
set yrange[0:0.25]

set xlabel '$t/ \tau_{\eta}$'
set ylabel '$\langle \zeta \rangle \qquad \textrm{var}(\zeta)$ '

set loadpath '../../../.config'
load 'spectral.pal'

set key

mhd = 0 

if(mhd == 0){

    load 'blues.pal'
    dir = "../data/sim_"
    name = "/line_evo_"
    file = ".dat"
    sim1 = "B05"
    file1 = dir.sim1.name.sim1.file

    set output "figures/hd_line_evo.png"
    plot file1 u 1:2 w l title sim1.': $\langle \zeta (t) \rangle$' ls 18,\
         file1 u 1:3 w l title sim1.': var($\zeta(t))$' ls 28,\
    }

    if(mhd == 1){

        load 'reds.pal'
        sim1 = "Q02"
        dir = "../data/sim_"
        name = "/line_evo_"
        file = ".dat"
        file1 = dir.sim1.name.sim1.file
        sim2 = "Q10"
        file2 = dir.sim2.name.sim2.file
        #sim3 = "Q02"
        #file3 = dir.sim3.name.sim3.file
        #sim4 = "Q10"
        #file4 = dir.sim4.name.sim4.file

        set output "figures/mhd_line_evo.tex"
        plot file1 u 1:2 w l title sim1.': $\langle \zeta (t) \rangle$' ls 18,\
             file1 u 1:3 w l title sim1.': var($\zeta(t))$' ls 28,\
             file2 u 1:2 w l title sim2.': $\langle \zeta (t) \rangle$' ls 17,\
             file2 u 1:3 w l title sim2.': var($\zeta(t))$' ls 27,\
             #file3 u 1:2 w l title sim3.': $\langle \zeta (t) \rangle$' ls 16,\
             #file3 u 1:3 w l title sim3.': var($\zeta(t))$' ls 26,\
             #file4 u 1:2 w l title sim4.': $\langle \zeta (t) \rangle$' ls 15,\
             #file4 u 1:3 w l title sim4.': var($\zeta(t))$' ls 25
}

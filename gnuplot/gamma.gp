reset

set terminal epslatex color size 6,4.5 

unset grid

set xrange[0:100]
set yrange[:]

set xlabel '$t/ \tau_{\eta}$'
set ylabel '$\langle \Gamma \rangle$'

set loadpath '../../../.config'
load 'spectral.pal'

mhd = 1 

if(mhd == 0){
    set key below

    load 'blues.pal'
    dir = "../data/sim_"
    name = "/gamma_"
    file = ".dat"
    sim1 = "C05"
    file1 = dir.sim1.name.sim1.file
    sim2 = "C07"
    file2 = dir.sim2.name.sim2.file

    set output "figures/hd_line_strain_evo.tex"

    plot file1 u 1:2 w l title sim1.': $\langle \Gamma(\va{T}_1,\va{l}) \rangle$' ls 18,\
         file1 u 1:3 w l title sim1.': $\langle \Gamma(\va{T}_2,\va{l}) \rangle$' ls 28,\
         file1 u 1:4 w l title sim1.': $\langle \Gamma(\va{T}_3,\va{l}) \rangle$' ls 38,\
         file1 u 1:5 w l title sim1.': $\langle \Gamma(\va{\omega},\va{l}) \rangle$' ls 48,\
         file2 u 1:2 w l title sim2.': $\langle \Gamma(\va{T}_1,\va{l}) \rangle$' ls 17,\
         file2 u 1:3 w l title sim2.': $\langle \Gamma(\va{T}_2,\va{l}) \rangle$' ls 27,\
         file2 u 1:4 w l title sim2.': $\langle \Gamma(\va{T}_3,\va{l}) \rangle$' ls 37,\
         file2 u 1:5 w l title sim2.': $\langle \Gamma(\va{\omega},\va{l}) \rangle$' ls 47,\


    set output "figures/hd_surf_strain_evo.tex"

    plot file1 u 1:12 w l title sim1.': $\langle \Gamma(\va{T}_1,\va{l}) \rangle$' ls 18,\
         file1 u 1:13 w l title sim1.': $\langle \Gamma(\va{T}_2,\va{l}) \rangle$' ls 28,\
         file1 u 1:14 w l title sim1.': $\langle \Gamma(\va{T}_3,\va{l}) \rangle$' ls 38,\
         file1 u 1:15 w l title sim1.': $\langle \Gamma(\va{\omega},\va{l}) \rangle$' ls 48,\
         file2 u 1:12 w l title sim2.': $\langle \Gamma(\va{T}_1,\va{l}) \rangle$' ls 17,\
         file2 u 1:13 w l title sim2.': $\langle \Gamma(\va{T}_2,\va{l}) \rangle$' ls 27,\
         file2 u 1:14 w l title sim2.': $\langle \Gamma(\va{T}_3,\va{l}) \rangle$' ls 37,\
         file2 u 1:15 w l title sim2.': $\langle \Gamma(\va{\omega},\va{l}) \rangle$' ls 47,\

}

if(mhd == 1){
    set key below

    load 'reds.pal'
    dir = "../data/sim_"
    name = "/gamma_"
    file = ".dat"
    sim1 = "Q10"
    file1 = dir.sim1.name.sim1.file
    sim2 = "Q02"
    file2 = dir.sim2.name.sim2.file

    set output "figures/mhd_line_strain_evo.tex"

    plot file1 u 1:2 w l title sim1.': $\langle \Gamma(\va{T}_1,\va{l}) \rangle$' ls 17,\
         file1 u 1:3 w l title sim1.': $\langle \Gamma(\va{T}_2,\va{l}) \rangle$' ls 27,\
         file1 u 1:4 w l title sim1.': $\langle \Gamma(\va{T}_3,\va{l}) \rangle$' ls 37,\
         file1 u 1:5 w l title sim1.': $\langle \Gamma(\va{\omega},\va{l}) \rangle$' ls 47,\
         file1 u 1:6 w l title sim1.': $\langle \Gamma(\va{b},\va{l}) \rangle$' ls 57,\
         file2 u 1:2 w l title sim2.': $\langle \Gamma(\va{T}_1,\va{l}) \rangle$' ls 18,\
         file2 u 1:3 w l title sim2.': $\langle \Gamma(\va{T}_2,\va{l}) \rangle$' ls 28,\
         file2 u 1:4 w l title sim2.': $\langle \Gamma(\va{T}_3,\va{l}) \rangle$' ls 38,\
         file2 u 1:5 w l title sim2.': $\langle \Gamma(\va{\omega},\va{l}) \rangle$' ls 48,\
         file2 u 1:6 w l title sim2.': $\langle \Gamma(\va{b},\va{l}) \rangle$' ls 58


    set output "figures/mhd_surf_strain_evo.tex"

    plot file1 u 1:12 w l title sim1.': $\langle \Gamma(\va{T}_1,\va{l}) \rangle$' ls 17,\
         file1 u 1:13 w l title sim1.': $\langle \Gamma(\va{T}_2,\va{l}) \rangle$' ls 27,\
         file1 u 1:14 w l title sim1.': $\langle \Gamma(\va{T}_3,\va{l}) \rangle$' ls 37,\
         file1 u 1:15 w l title sim1.': $\langle \Gamma(\va{\omega},\va{l}) \rangle$' ls 47,\
         file1 u 1:16 w l title sim1.': $\langle \Gamma(\va{b},\va{l}) \rangle$' ls 57,\
         file2 u 1:12 w l title sim2.': $\langle \Gamma(\va{T}_1,\va{l}) \rangle$' ls 18,\
         file2 u 1:13 w l title sim2.': $\langle \Gamma(\va{T}_2,\va{l}) \rangle$' ls 28,\
         file2 u 1:14 w l title sim2.': $\langle \Gamma(\va{T}_3,\va{l}) \rangle$' ls 38,\
         file2 u 1:15 w l title sim2.': $\langle \Gamma(\va{\omega},\va{l}) \rangle$' ls 48,\
         file2 u 1:16 w l title sim2.': $\langle \Gamma(\va{b},\va{l}) \rangle$' ls 58

}

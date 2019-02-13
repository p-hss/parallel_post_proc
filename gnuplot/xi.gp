reset

set terminal epslatex color size 6,4 

set loadpath '../../../.config'
load 'spectral.pal'

set yrange[0:0.25]
set xrange[0:100]

set xlabel '$t/ \tau_{\eta}$'
set ylabel '$\langle \xi \rangle \quad \textrm{var}(\xi)$ '

set key
unset grid

mhd = 0 

if(mhd == 0){

    load 'blues.pal'
    sim1 = "C05"
    dir = "../data/sim_"
    name = "/surf_evo_"
    file = ".dat"
    file1 = dir.sim1.name.sim1.file
    sim2 = "C07"
    file2 = dir.sim2.name.sim2.file
    sim3 = "C05"
    file3 = dir.sim3.name.sim3.file
    sim4 = "C07"
    file4 = dir.sim4.name.sim4.file

    set output "figures/hd_surf_evo.tex"

    plot file1 u 1:2 w l title sim1.': $\langle \xi (t) \rangle$' ls 18, \
         file1 u 1:3 w l title sim1.': $\textrm{var}(\xi (t))$' ls 28,\
         file2 u 1:2 w l title sim2.': $\langle \xi (t) \rangle$' ls 17,\
         file2 u 1:3 w l title sim2.': var($\xi(t))$' ls 27,\
         #file3 u 1:2 w l title sim3.': $\langle \xi (t) \rangle$' ls 16,\
         file3 u 1:3 w l title sim3.': var($\xi(t))$' ls 26,\
         file4 u 1:2 w l title sim4.': $\langle \xi (t) \rangle$' ls 16,\
         file4 u 1:3 w l title sim4.': var($\xi(t))$' ls 26
}

if(mhd == 1){

    load 'reds.pal'
    sim1 = "Q02"
    dir = "../data/sim_"
    name = "/surf_evo_"
    file = ".dat"
    file1 = dir.sim1.name.sim1.file
    sim2 = "Q10"
    file2 = dir.sim2.name.sim2.file
    sim3 = "Q02"
    file3 = dir.sim3.name.sim3.file
    sim4 = "Q10"
    file4 = dir.sim4.name.sim4.file

    set output "figures/mhd_surf_evo.tex"
    plot file1 u 1:2 w l  title sim1.': $\langle \xi (t) \rangle$' ls 18, \
         file1 u 1:3 w l  title sim1.': $\textrm{var}(\xi (t))$' ls 28,\
         file2 u 1:2 w l  title sim2.': $\langle \xi (t) \rangle$' ls 17, \
         file2 u 1:3 w l  title sim2.': $\textrm{var}(\xi (t))$' ls 27,\
         #file3 u 1:2 w l  title sim3.': $\langle \xi (t) \rangle$' ls 16, \
         file3 u 1:3 w l  title sim3.': $\textrm{var}(\xi (t))$' ls 26,\
         file4 u 1:2 w l  title sim4.': $\langle \xi (t) \rangle$' ls 15, \
         file4 u 1:3 w l  title sim4.': $\textrm{var}(\xi (t))$' ls 25
  }

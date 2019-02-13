reset

set terminal epslatex color size 6,4 

set xrange[:]
set yrange[0:]

set xlabel '$t/ \tau_{\eta}$'
set ylabel '$\langle |\va{v}_{\parallel}| \rangle \qquad \langle | \va{v}_{\perp} | \rangle$ '

set loadpath '../../../.config'
load 'spectral.pal'

set key bottom right

mhd = 1 

if(mhd == 0){

    load 'blues.pal'
    sim = "b05"
    dir = "../data/sim_"
    load_path = dir.sim
    name = "/mean_vel_"
    file = ".dat"
    file1 = dir.sim.name.sim.file

    set output "figures/hd_mean_vel.tex"
    plot file1 u 1:2 w l title\
         'B5: $\langle | \va{v}_{\parallel}| \rangle$' ls 18,\
         file1 u 1:3 w l title\
         'B5: $\langle | \va{v}_{\perp}| \rangle$' ls 28 
}

if(mhd == 1){

    load 'reds.pal'
    sim = "m02"
    dir = "../data/sim_"
    load_path = dir.sim
    name = "/mean_vel_"
    file = ".txt"
    file1 = dir.sim.name.sim.file

    sim = "m10"
    file = ".dat"
    file2 = dir.sim.name.sim.file

    set output "figures/mhd_mean_vel.tex"
    plot file1 u 1:2 w l title 'M2: $\langle | \va{v}_{\parallel}| \rangle$' ls 18,\
         file1 u 1:3 w l title 'M2: $\langle | \va{v}_{\perp}| \rangle$' ls 28,\
         file2 u 1:2 w l title 'M10: $\langle | \va{v}_{\parallel}| \rangle$' ls 17,\
         file2 u 1:3 w l title 'M10: $\langle | \va{v}_{\perp}| \rangle$' ls 27,\
}



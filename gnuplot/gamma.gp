reset

set terminal epslatex color size 6,4.5 

unset grid

set xrange[0:30]
set yrange[0.6:1.3]

set xlabel '$t/ \tau_{\eta}$'
set ylabel '$\langle$angle$\rangle$'

set loadpath '../../gnuplot_palettes/'
load 'spectral_poster.pal'
#load 'contrast.pal'

#set key below Left revers
set key below

sim = "M04 H00"
names = "MHD HD"

s1="100 1000"
s2="200 2000"
s3="300 3000"
s4="400 4000"
s5="500 5000"

files_1 = "../data/sim_"
files_2 = "/gamma_"
files_3 = ".dat"

set output "figures/line_strain_angle_evo.tex"
plot for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:2 w l \
     title word(names,i).': $\langle \sphericalangle(\va{T}_1,\va{l}) \rangle$' ls word(s1,i),\
     for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:3 w l \
     title word(names,i).': $\langle \sphericalangle(\va{T}_2,\va{l}) \rangle$' ls word(s2,i),\
     for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:4 w l \
     title word(names,i).': $\langle \sphericalangle(\va{T}_3,\va{l}) \rangle$' ls word(s3,i),\
     for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:5 w l \
     title word(names,i).': $\langle \sphericalangle(\va{\omega},\va{l}) \rangle$' ls word(s4,i),\
     files_1.word(sim,1).files_2.word(sim,1).files_3 u 1:6 w l \
     title word(names,1).': $\langle \sphericalangle(\va{b},\va{l}) \rangle$' ls word(s5,1)

#set yrange[0.6:1.5]
#
#files_1 = "../data/sim_"
#files_2 = "/gamma_"
#files_3 = ".dat"
#
#set output "figures/surface_strain_angle_evo.tex"
#plot for [i=1:words(sim)]\
#     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:12 w l \
#     title word(names,i).' $\langle \sphericalangle(\va{T}_1,\va{A}) \rangle$' ls i+9,\
#     for [i=1:words(sim)]\
#     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:13 w l \
#     title word(names,i).' $\langle \sphericalangle(\va{T}_2,\va{A}) \rangle$' ls i+19,\
#     for [i=1:words(sim)]\
#     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:14 w l \
#     title word(names,i).' $\langle \sphericalangle(\va{T}_3,\va{A}) \rangle$' ls i+29,\
#     for [i=1:words(sim)]\
#     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:15 w l \
#     title word(names,i).' $\langle \sphericalangle(\va{\omega},\va{A}) \rangle$' ls i+39,\
#     files_1.word(sim,1).files_2.word(sim,1).files_3 u 1:16 w l \
#     title word(names,1).' $\langle \sphericalangle(\va{b},\va{A}) \rangle$' ls 50
#
#mhd = 1 
#
#if(mhd == 0){
#    set key below
#
#    load 'blues.pal'
#    dir = "../data/sim_"
#    name = "/gamma_"
#    file = ".dat"
#    sim1 = "C05"
#    file1 = dir.sim1.name.sim1.file
#    sim2 = "C07"
#    file2 = dir.sim2.name.sim2.file
#
#    set output "figures/line_strain_angle_evo.tex"
#
#    plot file1 u 1:2 w l title sim1.': $\langle \Gamma(\va{T}_1,\va{l}) \rangle$' ls 18,\
#         file1 u 1:3 w l title sim1.': $\langle \Gamma(\va{T}_2,\va{l}) \rangle$' ls 28,\
#         file1 u 1:4 w l title sim1.': $\langle \Gamma(\va{T}_3,\va{l}) \rangle$' ls 38,\
#         file1 u 1:5 w l title sim1.': $\langle \Gamma(\va{\omega},\va{l}) \rangle$' ls 48,\
#         file2 u 1:2 w l title sim2.': $\langle \Gamma(\va{T}_1,\va{l}) \rangle$' ls 17,\
#         file2 u 1:3 w l title sim2.': $\langle \Gamma(\va{T}_2,\va{l}) \rangle$' ls 27,\
#         file2 u 1:4 w l title sim2.': $\langle \Gamma(\va{T}_3,\va{l}) \rangle$' ls 37,\
#         file2 u 1:5 w l title sim2.': $\langle \Gamma(\va{\omega},\va{l}) \rangle$' ls 47,\
#
#    set output "figures/hd_surf_strain_evo.tex"
#
#    plot file1 u 1:12 w l title sim1.': $\langle \Gamma(\va{T}_1,\va{l}) \rangle$' ls 18,\
#         file1 u 1:13 w l title sim1.': $\langle \Gamma(\va{T}_2,\va{l}) \rangle$' ls 28,\
#         file1 u 1:14 w l title sim1.': $\langle \Gamma(\va{T}_3,\va{l}) \rangle$' ls 38,\
#         file1 u 1:15 w l title sim1.': $\langle \Gamma(\va{\omega},\va{l}) \rangle$' ls 48,\
#         file2 u 1:12 w l title sim2.': $\langle \Gamma(\va{T}_1,\va{l}) \rangle$' ls 17,\
#         file2 u 1:13 w l title sim2.': $\langle \Gamma(\va{T}_2,\va{l}) \rangle$' ls 27,\
#         file2 u 1:14 w l title sim2.': $\langle \Gamma(\va{T}_3,\va{l}) \rangle$' ls 37,\
#         file2 u 1:15 w l title sim2.': $\langle \Gamma(\va{\omega},\va{l}) \rangle$' ls 47,\
#
#}
#
#if(mhd == 1){
#    set key below
#
#    load 'reds.pal'
#    dir = "../data/sim_"
#    name = "/gamma_"
#    file = ".dat"
#    sim1 = "Q10"
#    file1 = dir.sim1.name.sim1.file
#    sim2 = "Q02"
#    file2 = dir.sim2.name.sim2.file
#
#    set output "figures/mhd_line_strain_evo.tex"
#
#    plot file1 u 1:2 w l title sim1.': $\langle \Gamma(\va{T}_1,\va{l}) \rangle$' ls 17,\
#         file1 u 1:3 w l title sim1.': $\langle \Gamma(\va{T}_2,\va{l}) \rangle$' ls 27,\
#         file1 u 1:4 w l title sim1.': $\langle \Gamma(\va{T}_3,\va{l}) \rangle$' ls 37,\
#         file1 u 1:5 w l title sim1.': $\langle \Gamma(\va{\omega},\va{l}) \rangle$' ls 47,\
#         file1 u 1:6 w l title sim1.': $\langle \Gamma(\va{b},\va{l}) \rangle$' ls 57,\
#         file2 u 1:2 w l title sim2.': $\langle \Gamma(\va{T}_1,\va{l}) \rangle$' ls 18,\
#         file2 u 1:3 w l title sim2.': $\langle \Gamma(\va{T}_2,\va{l}) \rangle$' ls 28,\
#         file2 u 1:4 w l title sim2.': $\langle \Gamma(\va{T}_3,\va{l}) \rangle$' ls 38,\
#         file2 u 1:5 w l title sim2.': $\langle \Gamma(\va{\omega},\va{l}) \rangle$' ls 48,\
#         file2 u 1:6 w l title sim2.': $\langle \Gamma(\va{b},\va{l}) \rangle$' ls 58
#
#
#    set output "figures/mhd_surf_strain_evo.tex"
#
#    plot file1 u 1:12 w l title sim1.': $\langle \Gamma(\va{T}_1,\va{l}) \rangle$' ls 17,\
#         file1 u 1:13 w l title sim1.': $\langle \Gamma(\va{T}_2,\va{l}) \rangle$' ls 27,\
#         file1 u 1:14 w l title sim1.': $\langle \Gamma(\va{T}_3,\va{l}) \rangle$' ls 37,\
#         file1 u 1:15 w l title sim1.': $\langle \Gamma(\va{\omega},\va{l}) \rangle$' ls 47,\
#         file1 u 1:16 w l title sim1.': $\langle \Gamma(\va{b},\va{l}) \rangle$' ls 57,\
#         file2 u 1:12 w l title sim2.': $\langle \Gamma(\va{T}_1,\va{l}) \rangle$' ls 18,\
#         file2 u 1:13 w l title sim2.': $\langle \Gamma(\va{T}_2,\va{l}) \rangle$' ls 28,\
#         file2 u 1:14 w l title sim2.': $\langle \Gamma(\va{T}_3,\va{l}) \rangle$' ls 38,\
#         file2 u 1:15 w l title sim2.': $\langle \Gamma(\va{\omega},\va{l}) \rangle$' ls 48,\
#         file2 u 1:16 w l title sim2.': $\langle \Gamma(\va{b},\va{l}) \rangle$' ls 58
#
#}
#
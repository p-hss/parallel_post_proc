reset

set terminal epslatex color size 6,4 

set xrange[0:1]
set yrange[0:0.11]

set xlabel '$t/ \tau_{\eta}$'
set ylabel '$\langle \zeta \rangle \qquad \textrm{var}(\zeta)$ '

set loadpath '../../../.config'
load 'spectral.pal'

set key


#load 'blues.pal'
#dir = "../data/sim_"
#name = "/line_evo_"
#file = ".dat"
sim1 = "X02"
#file1 = dir.sim1.name.sim1.file
path_to_thesis = '/home/gamling/p.hess/Masterarbeit/Thesis/'

save path_to_thesis.'tables/tab_results_'.sim1.'.txt' file1x 
file1y = path_to_thesis.'tables/tab_mhd_sim_charact_'.sim1.'_1.txt'

print(file1x)

#set output "figures/mhd_cross_helicity_scaling.tex"
#plot "< cat -n file1x file1y" u 5:15 w p title sim1.': $\langle \zeta(\sigma_c) \rangle$' ls 18

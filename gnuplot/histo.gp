reset

set terminal epslatex color size 6,4 

set xrange[-1:1.3]
set yrange[:]

set loadpath '../../../.config'
load 'spectral.pal'

set xlabel '$\zeta$'
set ylabel 'Probability density function'
set key right top 

sim="R01 R02 R03 R04 R05 R06 R07 R08"

time="20"
files_1 = "../data/sim_"
files_2 = "/line_histo_"
files_3 = "_tkolmo_".time.".dat"

set output "figures/histograms/mhd_zeta_histo_t".time.".tex"
plot for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:2 w l \
     title word(sim,i).': $H(\zeta,\tau_{\eta}='.time.')$'  ls i+9

set xlabel '$\xi$'
set output "figures/histograms/mhd_xi_histo_t".time.".tex"
plot for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:2 w l \
     title word(sim,i).': $H(\xi,\tau_{\eta}='.time.'20)$'  ls i+9

reset

set terminal epslatex color size 6,4 

set xrange[0:80]

set xlabel '$t/ \tau_{\eta}$'
set key below

set loadpath '../../../.config'
load 'spectral.pal'

#sim = "C05 C01"
sim = "R01 R02 R03 R04 R05 R06 R07 R08"

files_1 = "../data/sim_"
files_2 = "/line_evo_"
files_3 = ".dat"

set yrange[-2*10**(-8):2*10**(-8)]
set ylabel '$\textrm{Skew}(\zeta)$'
set ytics 0.5*10**(-8)
set output "figures/line_higher_stats_skew_zeta.tex"
plot for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:5 w l \
     title word(sim,i).': $\textrm{Skew}(\zeta)$'  ls i+9

set ytics 0.5
set yrange[2:10]
set ylabel '$\textrm{Kurt}(\zeta)$'
set output "figures/line_higher_stats_kurt_zeta.tex"
plot for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:6 w l \
     title word(sim,i).': $\textrm{Kurt}(\zeta)$'  ls i+9

files_1 = "../data/sim_"
files_2 = "/surf_evo_"
files_3 = ".dat"

set yrange[-2*10**(-8):2*10**(-8)]
set ylabel '$\textrm{Skew}(\xi)$'
set ytics 0.5*10**(-8)
set output "figures/line_higher_stats_skew_xi.tex"
plot for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:5 w l \
     title word(sim,i).': $\textrm{Skew}(\xi)$'  ls i+9

set ytics 0.5
set yrange[2:10]
set ylabel '$\textrm{Kurt}(\xi)$'
set output "figures/line_higher_stats_kurt_xi.tex"
plot for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:6 w l \
     title word(sim,i).': $\textrm{Kurt}(\xi)$'  ls i+9

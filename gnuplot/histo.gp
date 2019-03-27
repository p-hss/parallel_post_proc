reset

set terminal epslatex color size 6,4 

set xrange[-1:1]
set yrange[0.001:2.5]

set loadpath '../../gnuplot_palettes/'
load 'spectral_poster.pal'

set ylabel 'Probability density function'
set key below 

sim="R01 R02 R03 R04 R05 R06"

angle="".sprintf("%1.2f", cos(0.1))." ".\
       sprintf("%1.2f", cos(0.2))." ".\
       sprintf("%1.2f", cos(0.3))." ".\
       sprintf("%1.2f", cos(0.4))." ".\
       sprintf("%1.2f", cos(0.5))." ".\
       sprintf("%1.2f", cos(0.6))." "

time="20"
files_1 = "../data/sim_"
files_2 = "/line_histo_"
files_3 = "_tkolmo_".time.".dat"

#set xtics 10
set logscale y
#set logscale x
#set format y "$10^{%L}$"


G1(x) = 1./(sigma1*sqrt(2*pi)) * exp( -(x-mu1)**2 / (2*sigma1**2) )
G2(x) = 1./(sigma2*sqrt(2*pi)) * exp( -(x-mu2)**2 / (2*sigma2**2) )
#
zeta="0.115 0.110 0.106 0.098 0.101 0.0998"

fit G1(x) files_1.word(sim,1).files_2.word(sim,1).files_3 u ($1-word(zeta,1)):2 via sigma1, mu1
fit G2(x) files_1.word(sim,6).files_2.word(sim,6).files_3 u ($1-word(zeta,6)):2 via sigma2, mu2

set xlabel '$\zeta$'
set output "figures/histograms/mhd_zeta_histo_t".time.".tex"
plot for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u ($1-word(zeta,i)):2 w l \
     title ' $\sigma_{c, \textsf{f}}=$'.word(angle,i) ls i+9,\
     G1(x) ls 200  title "Gauss",\
     G2(x) ls 200  notitle 
     #title word(sim,i).': $H(\zeta,\tau_{\eta}='.time.')$'  ls i+9

set xlabel '$\xi$'
set output "figures/histograms/mhd_xi_histo_t".time.".tex"
plot for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:2 w l \
     title word(sim,i).': $H(\xi,\tau_{\eta}='.time.'20)$'  ls i+9,\


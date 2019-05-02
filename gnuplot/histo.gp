reset

set terminal epslatex color size 6,4 

set xrange[-10:10]
set yrange[0.001:5.0]

set loadpath '../../gnuplot_palettes/'
load 'spectral_poster.pal'

set ylabel 'Probability density function'
set key below Left reverse
set bmargin 7

#sim="Z01 Z02 Z03 Z04 Z05 Z06"
sim="Z41 Z42 Z43 Z44 Z45 Z46"

angle="".sprintf("%1.2f", 0.98)." ".\
         sprintf("%1.2f", 0.94)." ".\
         sprintf("%1.2f", 0.84)." ".\
         sprintf("%1.2f", 0.71)." ".\
         sprintf("%1.1f", 0.5)." ".\
         sprintf("%1.0f", 0)." "

time="20"
files_1 = "../data/sim_"
files_2 = "/line_histo_"
files_3 = "_tkolmo_".time.".dat"

set xtics 2
#set logscale y
#set logscale x
#set format y "$10^{%L}$"


G1(x) = 1./(sigma1*sqrt(2*pi)) * exp( -(x-mu1)**2 / (2*sigma1**2) )
G2(x) = 1./(sigma2*sqrt(2*pi)) * exp( -(x-mu2)**2 / (2*sigma2**2) )
G(x) = 1./(sqrt(2*pi)) * exp( -(x)**2/2 )
#
#zeta="0.118 0.116 0.109 0.100 0.098 0.099"
zeta="0.08 0.075 0.066 0.062 0.062 0.061"
sigma="0.205 0.189 0.169 0.161 0.157 0.152"

#fit G(x) files_1.word(sim,6).files_2.word(sim,6).files_3 u ($1-word(zeta,6))/word(sigma,6):2 via a
#fit G1(x) files_1.word(sim,1).files_2.word(sim,1).files_3 u ($1-word(zeta,1))/word(sigma,1):2 via sigma1, mu1
#fit G2(x) files_1.word(sim,6).files_2.word(sim,6).files_3 u ($1-word(zeta,6)):2 via sigma2, mu2

set xlabel '$(\zeta-\langle \zeta \rangle)/\sigma$'
#set xlabel offset 0,0.5

set output "figures/histograms/mhd_zeta_Z4_2_histo_t".time."_256.tex"
plot for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u ($1-word(zeta,i))/word(sigma,i):2 w l \
     title ' $\sigma_{c}=$ '.word(angle,i) ls i+9,\
     G(x) ls 200  title "Gauss"
     #G2(x) ls 200  notitle 
     #title word(sim,i).': $H(\zeta,\tau_{\eta}='.time.')$'  ls i+9

set xlabel '$\xi$'
set output "figures/histograms/mhd_xi_Z4_histo_t".time.".tex"
plot for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:2 w l \
     title word(sim,i).': $H(\xi,\tau_{\eta}='.time.'20)$'  ls i+9,

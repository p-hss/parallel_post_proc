reset

set terminal epslatex color size 6,3.5 
#set terminal pdf

set xrange[-5:5]
set yrange[0.0001:]

set loadpath '../../gnuplot_palettes/'
load 'spectral_poster.pal'
#load 'contrast.pal'

set ylabel 'Probability density function'
#set key below Left reverse
set key

sim="M04 H00"
names="MHD HD"
zeta="0.069 0.112"
sigma="0.095 0.123" 

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

s1="100 1000"
s2="200 2000"

set xtics 2
set logscale y
#set format y "$10^{%L}$"

G1(x) = 1./(sigma1*sqrt(2*pi)) * exp( -(x-mu1)**2 / (2*sigma1**2) )
G2(x) = 1./(sigma2*sqrt(2*pi)) * exp( -(x-mu2)**2 / (2*sigma2**2) )
G(x) = 1./(sqrt(2*pi)) * exp( -(x)**2/2 )

#zeta="0.118 0.116 0.109 0.100 0.098 0.099" #Z0
#zeta="0.08 0.075 0.066 0.062 0.062 0.061" #Z4
#sigma="0.205 0.189 0.169 0.161 0.157 0.152" #Z4
#zeta="0.1217"
#sigma="0.182"

#fit G(x) files_1.word(sim,6).files_2.word(sim,6).files_3 u ($1-word(zeta,6))/word(sigma,6):2 via a
fit [-1:1] G1(x) files_1.word(sim,1).files_2.word(sim,1).files_3 u 1:2 via sigma1, mu1
fit [-1:1] G2(x) files_1.word(sim,2).files_2.word(sim,2).files_3 u 1:2 via sigma2, mu2

set xlabel '$(\zeta-\langle \zeta \rangle)/\sigma$'
#set xlabel offset 0,0.5

set output "figures/histograms/mhd_zeta_histo_new_t".time."_256.tex"
plot for [i=1:2]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u ($1-word(zeta,i))/word(sigma,i):($2*word(sigma,i)) w l \
     title word(names,i).'' ls word(s1,i),\
     G1(x*word(sigma,1)+word(zeta,1))*word(sigma,1) ls 200  title "Gauss fit"
     #G2(x*word(sigma,2)+word(zeta,2))*word(sigma,2) ls 200  notitle 
     #title word(sim,i).': $H(\zeta,\tau_{\eta}='.time.')$'  ls i+9

#set output "figures/histograms/mhd_zeta_histo_t".time."_256.tex"
#plot for [i=1:words(sim)]\
#     files_1.word(sim,i).files_2.word(sim,i).files_3 u ($1-word(zeta,i))/word(sigma,i):($2*word(sigma,i)) w l \
#     title word(names,i).'' ls i+9,\
#     #G(x) ls 200  title "Gauss",\
#     G2(x) ls 200  notitle 
#     #title word(sim,i).': $H(\zeta,\tau_{\eta}='.time.')$'  ls i+9

set xlabel '$\xi$'
set output "figures/histograms/mhd_xi_histo_t".time.".tex"
plot for [i=1:words(sim)]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u 1:2 w l \
     title word(sim,i).': $H(\xi,\tau_{\eta}='.time.'20)$'  ls i+9,


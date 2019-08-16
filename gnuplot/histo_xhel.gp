reset

set terminal epslatex color size 6,3.5 

set xrange[-5:5]
set yrange[0.0001:]

set loadpath '../../gnuplot_palettes/'
load 'spectral_poster.pal'
#load 'contrast.pal'

set ylabel 'PDF'
#set key below Left reverse
set key

sim="X40 X46"
names="XHEL1 XHEL6"
zeta="0.115 0.042"
sigma="0.074 0.024" 

zeta="0.116 0.045"
sigma="0.272 0.158" 

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

#fit G(x) files_1.word(sim,6).files_2.word(sim,6).files_3 u ($1-word(zeta,6))/word(sigma,6):2 via a
fit [-1:1] G1(x) files_1.word(sim,1).files_2.word(sim,1).files_3 u 1:2 via sigma1, mu1
fit [-1:1] G2(x) files_1.word(sim,2).files_2.word(sim,2).files_3 u 1:2 via sigma2, mu2

set xlabel '$(\zeta-\langle \zeta \rangle)/\sigma$'
#set xlabel offset 0,0.5

set output "figures/histograms/xhel_zeta_histo_new_t".time."_256.tex"
plot for [i=1:2]\
     files_1.word(sim,i).files_2.word(sim,i).files_3 u ($1-word(zeta,i))/word(sigma,i):($2*word(sigma,i)) w l \
     title word(names,i).'' ls word(s1,i),\
     G1(x*word(sigma,1)+word(zeta,1))*word(sigma,1) ls 200  title "Gauss fit"
     #G2(x*word(sigma,2)+word(zeta,2))*word(sigma,2) ls 200  notitle 
     #title word(sim,i).': $H(\zeta,\tau_{\eta}='.time.')$'  ls i+9

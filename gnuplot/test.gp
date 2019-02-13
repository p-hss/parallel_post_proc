reset

set terminal epslatex color size 6,4 


set xrange[0:35]
set yrange[:]

set xlabel '$t/ \tau_{\eta}$'

set loadpath '../../../.config'
load 'spectral.pal'
set palette negative

unset arrow
unset key

set loadpath '../data'
filename = "line_test.dat"

set ylabel '$\ln(l)$'
set output "figures/single_le.tex"
plot for [i=33:52] filename using 1:i:i-31 w l palette lt 1 lw 4

unset colorbox
set ylabel '$\zeta$'
set output "figures/single_zeta.tex"
plot for [i=3:22] filename using 1:i:(i+21) w l palette lt 1 lw 4


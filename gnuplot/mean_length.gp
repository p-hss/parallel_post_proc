reset

set terminal epslatex color size 6,4 

set xrange[0:]
set yrange[:]

set xlabel '$t/ \tau_{\eta}$'
#set logscale y

set loadpath '../../../.config'
load 'spectral.pal'

set key right bottom 
unset grid

sim = "b05"
sim_up = "B05"
dir = "../data/sim_"
set loadpath dir.sim
name = "/mean_length_evo_"
file = ".dat"
file1 = name.sim.file

line_fit(x) = zeta*x + b
surface_fit(x) = xi*x + d

fit [10:35] line_fit(x) file1 using 1:2:3 yerr via zeta,b 
fit [10:35] surface_fit(x) file1 using 1:4 via xi, d 

set ylabel '$\langle \ln(l)\rangle$'
set output "figures/mean_line_evo.tex"
plot file1 u 1:2:3 w errorbars title 'B5: $\langle \ln(l(t)) \rangle$' ls 40,\
     line_fit(x) title 'linear fit: $\zeta x + c $' ls 1 lw 5

set ylabel '$\langle \ln(A)\rangle$'
set output "figures/mean_surf_evo.tex"
plot file1 u 1:4:5 w errorbars title 'B5: $\langle \ln(A(t)) \rangle$' ls 40,\
     surface_fit(x) title 'linear fit: $\xi x + c $' ls 1 lw 5

#---------------------Output to file
fit_paramter_file = "../tables/tab_growth_rates_"
set print fit_paramter_file.sim.file 
lb = "\\"
print sprintf('%s & %2.3f & %1.3E & %2.3f & %1.3E', sim_up, zeta, zeta_err, xi,\
xi_err)  




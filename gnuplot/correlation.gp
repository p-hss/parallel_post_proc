reset

set terminal epslatex color size 6,4 

set key

set yrange[:]
set xrange[0:20]

set xlabel '$\tau_{\textrm{lag}}=t-t^{*}$'
set ylabel '$\textrm{cov}(\zeta(t)\zeta(t+\tau_{\textrm{lag}})) \quad \
\textrm{cov}(\xi(t)\xi(t+\tau_{\textrm{lag}}))$'
set zeroaxis

set loadpath '../../../.config'
load 'spectral.pal'

mhd = 1 

if(mhd == 0){
    load 'blues.pal'
    sim = "b05"
    dir = "../data/sim_"
    set loadpath dir.sim
    name1 = "/covariance_"
    file = ".dat"
    file1 = name1.sim.file

    set output "figures/hd_covariance.tex"

    plot file1 u (($1)*1.8/1.2):3 w l title 'B5: \
    $\textrm{cov}(\xi(t)\xi(t+\tau_{\textrm{lag}}))$' ls 18,\
         file1 u 1:2 w l title 'B5: \
    $\textrm{cov}(\zeta(t)\zeta(t+\tau_{\textrm{lag}}))$' ls 28
}

if(mhd == 1){

    load 'reds.pal'

    sim = "m02"
    dir = "../data/sim_"
    set loadpath dir.sim
    name1 = "/covariance_"
    file = ".dat"
    file1 = name1.sim.file

    set output "figures/mhd_covariance.tex"

    plot file1 u 1:3 w l title 'M2: \
    $\textrm{cov}(\xi(t)\xi(t+\tau_{\textrm{lag}}))$' ls 18,\
         file1 u 1:2 w l title 'M2: \
    $\textrm{cov}(\zeta(t)\zeta(t+\tau_{\textrm{lag}}))$' ls 28
}

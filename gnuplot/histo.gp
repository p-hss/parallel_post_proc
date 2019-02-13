reset

#set terminal epslatex color size 6,4 
set terminal png

set xrange[-1:1.3]
set yrange[:]

set loadpath '../../../.config'
load 'spectral.pal'

set xlabel '$\zeta$'
set ylabel 'Probability density function'

set key right top 

mhd = 0 

if(mhd == 0){
    load 'blues.pal'

    dir = "../data/sim_"
    name1 = "/line_histo_"
    name2 = "_tkolmo_20"
    file = ".dat"
    sim1 = "B05"
    file1 = dir.sim1.name1.sim1.name2.file

    set output "figures/histograms/hd_zeta_histo_t10.png"
    plot file1 u 1:2 w l  title sim1.': $H(\zeta,\tau_{\eta}=10)$' ls 18

    #name1 = "/surf_histo_"
    #file1 = dir.sim1.name1.sim1.name2.file
    #file2 = dir.sim2.name1.sim2.name2.file

    #set xlabel '$\xi$'
    #set output "figures/histograms/hd_xi_histo_t10.tex"
    #plot file1 u 1:2 w l  title sim1.': $H(\xi,\tau_{\eta}=10)$' ls 18,\
    #     file2 u 1:2 w l  title sim1.': $H(\xi,\tau_{\eta}=10)$' ls 17 

}

if(mhd == 1){

    load 'reds.pal'

    dir = "../data/sim_"
    name1 = "/line_histo_"
    name2 = "_tkolmo_20"
    file = ".dat"
    sim1 = "Q02"
    file1 = dir.sim1.name1.sim1.name2.file
    sim2 = "Q10"
    file2 = dir.sim2.name1.sim2.name2.file

    set output "figures/histograms/mhd_zeta_histo_t20.tex"
    plot file1 u 1:2 w l  title sim1.': $H(\zeta,\tau_{\eta}=20)$' ls 18,\
         file2 u 1:2 w l  title sim2.': $H(\zeta,\tau_{\eta}=20)$' ls 17

    name1 = "/surf_histo_"
    file1 = dir.sim1.name1.sim1.name2.file
    file2 = dir.sim2.name1.sim2.name2.file

    set xlabel '$\xi$'
    set output "figures/histograms/mhd_xi_histo_t20.tex"
    plot file1 u 1:2 w l  title sim1.': $H(\xi,\tau_{\eta}=20)$' ls 18,\
         file2 u 1:2 w l  title sim2.': $H(\xi,\tau_{\eta}=20)$' ls 17
}

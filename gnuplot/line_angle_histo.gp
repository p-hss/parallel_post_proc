reset

set terminal epslatex color size 6,4.8 

set xrange[:]
set yrange[:]

set loadpath '../../../.config'

set ylabel 'Histogram'


mhd = 1 

if(mhd == 0){
    
    set key left top 

    load 'blues.pal'
    dir = "../data/sim_"
    name1 = "/angle_histo_"
    name2 = "_tkolmo_20"

    sim1 = "C05"
    file = ".dat"
    file1 = dir.sim1.name1.sim1.name2.file
    sim2 = "C07"
    file = ".dat"
    file2 = dir.sim2.name1.sim2.name2.file

    set output "figures/histograms/hd_angle_histo_t20.tex"
    set xlabel '$|\cos(\Gamma(\va{T}_i,\va{l}))|$'

    plot file1 u 1:2 w l \
        title sim1.': $H(|\cos(\measuredangle(\va{T}_1,\va{l}))|)$' ls 18, \
        file1 u 1:3 w l \
        title sim1.': $H(|\cos(\measuredangle(\va{T}_2,\va{l}))|)$' ls 28, \
        file1 u 1:4 w l \
        title sim1.': $H(|\cos(\measuredangle(\va{T}_3,\va{l}))|)$' ls 38, \
        file1 u 1:5 w l\
        title sim1.': $H(|\cos(\measuredangle(\va{\omega},\va{l}))|)$' ls 48,\
        file2 u 1:2 w l \
        title sim2.': $H(|\cos(\measuredangle(\va{T}_1,\va{l}))|)$' ls 17, \
        file2 u 1:3 w l \
        title sim2.': $H(|\cos(\measuredangle(\va{T}_2,\va{l}))|)$' ls 27, \
        file2 u 1:4 w l \
        title sim2.': $H(|\cos(\measuredangle(\va{T}_3,\va{l}))|)$' ls 37, \
        file2 u 1:5 w l\
        title sim2.': $H(|\cos(\measuredangle(\va{\omega},\va{l}))|)$' ls 47
}

if(mhd == 1){
    set terminal epslatex color size 6,5 
    set key below

    load 'reds.pal'
    dir = "../data/sim_"
    name1 = "/angle_histo_"
    name2 = "_tkolmo_20"
    file = ".dat"

    sim1 = "Q02"
    file1 = dir.sim1.name1.sim1.name2.file
    file12 = dir.sim1."/mhd_angle_histo_".sim1.name2.file

    sim2 = "Q10"
    file2 = dir.sim2.name1.sim2.name2.file
    file22 = dir.sim2."/mhd_angle_histo_".sim2.name2.file 

    set output "figures/histograms/mhd_angle_histo_t20.tex"
    set xlabel '$|\cos(\Gamma(\va{T}_i,\va{l}))|$'

    plot file1 u 1:2 w l \
         title sim1.': $H(|\cos(\measuredangle(\va{T}_1,\va{l}))|)$' ls 18, \
         file1 u 1:3 w l \
         title sim1.': $H(|\cos(\measuredangle(\va{T}_2,\va{l}))|)$'ls 28,\
         file1 u 1:4 w l \
         title sim1.': $H(|\cos(\measuredangle(\va{T}_3,\va{l}))|)$'ls 38,\
         file1 u 1:5 w l \
         title sim1.': $H(|\cos(\measuredangle(\va{\omega},\va{l}))|)$' ls 48,\
         file12 u 1:2 w l \
         title sim1.': $H(|\cos(\measuredangle(\delta \va{b},\va{l}))|)$' ls 58,\
         file2 u 1:2 w l \
         title sim2.': $H(|\cos(\measuredangle(\va{T}_1,\va{l}))|)$' ls 17, \
         file2 u 1:3 w l \
         title sim2.': $H(|\cos(\measuredangle(\va{T}_2,\va{l}))|)$'ls 27,\
         file2 u 1:4 w l \
         title sim2.': $H(|\cos(\measuredangle(\va{T}_3,\va{l}))|)$'ls 37,\
         file2 u 1:5 w l \
         title sim2.': $H(|\cos(\measuredangle(\va{\omega},\va{l}))|)$' ls 47,\
         file22 u 1:2 w l \
         title sim2.': $H(|\cos(\measuredangle(\delta \va{b},\va{l}))|)$' ls 57
}

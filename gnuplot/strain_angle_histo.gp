reset

#set terminal epslatex color size 6,4 
set terminal png

set xrange[:]
set yrange[:]

set loadpath '../../../.config'

set ylabel 'Histogram'

set key top left

mhd = 0 

if(mhd == 0){

    load 'blues.pal'
    dir = "../data/sim_"
    name1 = "/angle_histo_"
    name2 = "_tkolmo_20"
    file = ".dat"

    sim1 = "B05"
    file1 = dir.sim1.name1.sim1.name2.file

    set xlabel '$|\cos(\Gamma(\va{T}_i,\va{\omega}))|$'
    #set output "figures/histograms/hd_strain_vorticity_angle_histo_t20.tex"
    set output "figures/histograms/hd_strain_vorticity_angle_histo_t20.png"

    plot file1 u 1:6 w l\
         title sim1.': $H(|\cos(\measuredangle(\va{T_1},\va{\omega}))|)$' ls 18, \
         file1 u 1:7 w l \
         title sim1.': $H(|\cos(\measuredangle(\va{T_2},\va{\omega}))|)$' ls 28, \
         file1 u 1:8 w l \
         title sim1.': $H(|\cos(\measuredangle(\va{T_3},\va{\omega}))|)$' ls 38
}

if(mhd == 1){
    load 'reds.pal'
    dir = "../data/sim_"
    name1 = "/mhd_angle_histo_"
    name2 = "_tkolmo_20"
    file = ".dat"


    sim1 = "Q10"
    file1 = dir.sim1.name1.sim1.name2.file

    sim2 = "Q02"
    file2 = dir.sim2.name1.sim2.name2.file

    set xlabel '$|\cos(\measuredangle(\va{b},\va{T}_i))|$'
    set output "figures/histograms/mhd_strain_magnetic_angle_histo_t20.tex"

    plot file1 u 1:3 w l \
         title sim1.': $H(|\cos(\measuredangle(\delta \va{b},\va{T_1}))|)$' ls 18, \
         file1 u 1:4 w l \
         title sim1.': $H(|\cos(\measuredangle(\delta \va{b},\va{T_2}))|)$' ls 28, \
         file1 u 1:5 w l \
         title sim1.': $H(|\cos(\measuredangle(\delta \va{b},\va{T_3}))|)$' ls 38, \
         file1 u 1:6 w l \
         title sim1.': $H(|\cos(\measuredangle(\delta \va{b},\va{\omega}))|)$' ls 48, \
         file2 u 1:3 w l \
         title sim2.': $H(|\cos(\measuredangle(\delta \va{b},\va{T_1}))|)$' ls 17, \
         file2 u 1:4 w l \
         title sim2.': $H(|\cos(\measuredangle(\delta \va{b},\va{T_2}))|)$' ls 27, \
         file2 u 1:5 w l \
         title sim2.': $H(|\cos(\measuredangle(\delta \va{b},\va{T_3}))|)$' ls 37, \
         file2 u 1:6 w l \
         title sim2.': $H(|\cos(\measuredangle(\delta \va{b},\va{\omega}))|)$' ls 47, \


    set xlabel '$|\cos(\Gamma(\va{T}_i,\va{\omega}))|$'
    set output "figures/histograms/mhd_strain_vorticity_angle_histo_t20.tex"

    plot file1 u 1:7 w l\
         title sim1.': $H(|\cos(\measuredangle(\va{T_1},\va{\omega}))|)$' ls 18, \
         file1 u 1:8 w l \
         title sim1.': $H(|\cos(\measuredangle(\va{T_2},\va{\omega}))|)$' ls 28, \
         file1 u 1:9 w l \
         title sim1.': $H(|\cos(\measuredangle(\va{T_3},\va{\omega}))|)$' ls 38, \
         file2 u 1:7 w l\
         title sim2.': $H(|\cos(\measuredangle(\va{T_1},\va{\omega}))|)$' ls 17, \
         file2 u 1:8 w l \
         title sim2.': $H(|\cos(\measuredangle(\va{T_2},\va{\omega}))|)$' ls 27, \
         file2 u 1:9 w l \
         title sim2.': $H(|\cos(\measuredangle(\va{T_3},\va{\omega}))|)$' ls 37, \
}

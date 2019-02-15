#!/bin/bash
#shell script tp run the code
make clean
make
#run the executable with 4 proccesses
mpirun -np 4 run < param/C05_param.dat 

#!/bin/bash
#shell script tp run the code
make clean
make
srun -n 2 ./run < param/Z01.dat 

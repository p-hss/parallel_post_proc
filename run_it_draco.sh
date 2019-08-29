#!/bin/bash
#shell script tp run the code
make clean
make
srun -n 2 ./run < param/B05.dat 

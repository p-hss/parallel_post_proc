#!/bin/bash

# plot growth rates
gnuplot zeta.gp
gnuplot xi.gp

#plot the time evolution of the angles 
gnuplot gamma.gp
gnuplot vol.gp

#eigenvalues of the cauchy green tensor
gnuplot cauchy_green.gp



# Parallelized code for computing material line statistics from MHD simulations

## Requirements
The `LAPACK` and `BLAS` libraries need to be installed and their paths set in the `make` file.

## How to use
Compile the code using the `make` file. The executable can then be run on a SLURM cluster, e.g. with
```
srun -n <N-nodes> ./run < param/<parameter-file-name>.dat 
```
where `<N-nodes>` is the number of nodes, e.g. 2, and `<parameter-file-name>.dat` is a configuration file in the `param/` directory. 
The input files from the MHD code are expected in  `../results/tmp/` and the output of the analysis is written to the `data/` directory. 
Further examples for running the code can be found in the `jobscripts` directory.

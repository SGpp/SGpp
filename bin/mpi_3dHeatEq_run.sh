#!/bin/sh

mpiexec -ppn 2 -genv OMP_NUM_THREADS 4 -genv I_MPI_PIN_DOMAIN omp:platform -n 4 ./HESolverMPI HeatEquation 3 11 0.0 3.0 1.0 smooth 1.0 0.1 ImEul 0.00001 400

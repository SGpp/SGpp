#!/bin/sh

mpiexec -ppn 2 -genv OMP_NUM_THREADS 4 -genv I_MPI_PIN_DOMAIN omp:platform -n 4 ./HESolverMPI PoissonEquation 4 10 0.0 3.0 smooth 0.00001 400

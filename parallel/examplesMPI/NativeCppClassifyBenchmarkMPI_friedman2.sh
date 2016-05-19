mpiexec -n 4 ./NativeCppClassifyBenchmarkMPI friedman2 friedman_4d_2000.arff 0 DP linearboundary 3 0.000001 250 0.0001 6 0.0 100 20 0.1 X86SIMD Bigdata
#possible values for mpi communication: Allreduce, Alltoallv, Async, Onesided, TrueAsync, Bigdata

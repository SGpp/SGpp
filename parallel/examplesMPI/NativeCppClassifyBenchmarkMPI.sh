mpiexec -n 28 ./NativeCppClassifyBenchmarkMPI chess_5d_2000.arff chess_train_5d_very_small.arff 0 DP linearboundary 3 0.000001 250 0.0001 6 0.0 100 20 0.1 X86SIMD Allreduce
#possible values for mpi communication: Allreduce, Alltoallv, Async, Onesided, TrueAsync, Bigdata

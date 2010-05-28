#!/bin/sh

OMP_NUM_THREADS=8
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 5 2 BStest5DMC.bound BStest5DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test5D_level_2_8CPU.out
./NativeCppBlackScholesSolver solveND 5 3 BStest5DMC.bound BStest5DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test5D_level_3_8CPU.out
./NativeCppBlackScholesSolver solveND 5 4 BStest5DMC.bound BStest5DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test5D_level_4_8CPU.out
./NativeCppBlackScholesSolver solveND 5 5 BStest5DMC.bound BStest5DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test5D_level_5_8CPU.out
./NativeCppBlackScholesSolver solveND 5 6 BStest5DMC.bound BStest5DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test5D_level_6_8CPU.out
./NativeCppBlackScholesSolver solveND 5 7 BStest5DMC.bound BStest5DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test5D_level_7_8CPU.out
OMP_NUM_THREADS=1
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 5 2 BStest5DMC.bound BStest5DMC.stoch 1.0 std_euro_call 0.00 1.0 1.0 CrNic 16000 0.00001 > test5D_level_2_1CPU.out
./NativeCppBlackScholesSolver solveND 5 3 BStest5DMC.bound BStest5DMC.stoch 1.0 std_euro_call 0.00 1.0 1.0 CrNic 16000 0.00001 > test5D_level_3_1CPU.out
./NativeCppBlackScholesSolver solveND 5 4 BStest5DMC.bound BStest5DMC.stoch 1.0 std_euro_call 0.00 1.0 1.0 CrNic 16000 0.00001 > test5D_level_4_1CPU.out
./NativeCppBlackScholesSolver solveND 5 5 BStest5DMC.bound BStest5DMC.stoch 1.0 std_euro_call 0.00 1.0 1.0 CrNic 16000 0.00001 > test5D_level_5_1CPU.out
./NativeCppBlackScholesSolver solveND 5 6 BStest5DMC.bound BStest5DMC.stoch 1.0 std_euro_call 0.00 1.0 1.0 CrNic 16000 0.00001 > test5D_level_6_1CPU.out
./NativeCppBlackScholesSolver solveND 5 7 BStest5DMC.bound BStest5DMC.stoch 1.0 std_euro_call 0.00 1.0 1.0 CrNic 16000 0.00001 > test5D_level_7_1CPU.out
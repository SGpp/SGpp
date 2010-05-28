#!/bin/sh

OMP_NUM_THREADS=8
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 4 2 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_2_8CPU.out
./NativeCppBlackScholesSolver solveND 4 3 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_3_8CPU.out
./NativeCppBlackScholesSolver solveND 4 4 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_4_8CPU.out
./NativeCppBlackScholesSolver solveND 4 5 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_5_8CPU.out
./NativeCppBlackScholesSolver solveND 4 6 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_6_8CPU.out
./NativeCppBlackScholesSolver solveND 4 7 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_7_8CPU.out
OMP_NUM_THREADS=4
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 4 2 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_2_4CPU.out
./NativeCppBlackScholesSolver solveND 4 3 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_3_4CPU.out
./NativeCppBlackScholesSolver solveND 4 4 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_4_4CPU.out
./NativeCppBlackScholesSolver solveND 4 5 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_5_4CPU.out
./NativeCppBlackScholesSolver solveND 4 6 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_6_4CPU.out
./NativeCppBlackScholesSolver solveND 4 7 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_7_4CPU.out
OMP_NUM_THREADS=2
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 4 2 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_2_2CPU.out
./NativeCppBlackScholesSolver solveND 4 3 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_3_2CPU.out
./NativeCppBlackScholesSolver solveND 4 4 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_4_2CPU.out
./NativeCppBlackScholesSolver solveND 4 5 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_5_2CPU.out
./NativeCppBlackScholesSolver solveND 4 6 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_6_2CPU.out
./NativeCppBlackScholesSolver solveND 4 7 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_7_2CPU.out
OMP_NUM_THREADS=1
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 4 2 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_2_1CPU.out
./NativeCppBlackScholesSolver solveND 4 3 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_3_1CPU.out
./NativeCppBlackScholesSolver solveND 4 4 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_4_1CPU.out
./NativeCppBlackScholesSolver solveND 4 5 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_5_1CPU.out
./NativeCppBlackScholesSolver solveND 4 6 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_6_1CPU.out
./NativeCppBlackScholesSolver solveND 4 7 BStest4DMC.bound BStest4DMC.stoch 1.0 std_euro_call 0.00 1.0 0.1 CrNic 16000 0.00001 > test4D_level_7_1CPU.out

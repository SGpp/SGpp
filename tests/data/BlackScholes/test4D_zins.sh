#!/bin/sh

OMP_NUM_THREADS=8
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 4 2 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_2_zins_8CPU.out
./NativeCppBlackScholesSolver solveND 4 3 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_3_zins_8CPU.out
./NativeCppBlackScholesSolver solveND 4 4 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_4_zins_8CPU.out
./NativeCppBlackScholesSolver solveND 4 5 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_5_zins_8CPU.out
./NativeCppBlackScholesSolver solveND 4 6 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_6_zins_8CPU.out
./NativeCppBlackScholesSolver solveND 4 7 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_7_zins_8CPU.out
OMP_NUM_THREADS=2
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 4 2 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_2_zins_4CPU.out
./NativeCppBlackScholesSolver solveND 4 3 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_3_zins_4CPU.out
./NativeCppBlackScholesSolver solveND 4 4 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_4_zins_4CPU.out
./NativeCppBlackScholesSolver solveND 4 5 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_5_zins_4CPU.out
./NativeCppBlackScholesSolver solveND 4 6 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_6_zins_4CPU.out
./NativeCppBlackScholesSolver solveND 4 7 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_7_zins_4CPU.out
OMP_NUM_THREADS=2
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 4 2 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_2_zins_2CPU.out
./NativeCppBlackScholesSolver solveND 4 3 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_3_zins_2CPU.out
./NativeCppBlackScholesSolver solveND 4 4 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_4_zins_2CPU.out
./NativeCppBlackScholesSolver solveND 4 5 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_5_zins_2CPU.out
./NativeCppBlackScholesSolver solveND 4 6 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_6_zins_2CPU.out
./NativeCppBlackScholesSolver solveND 4 7 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_7_zins_2CPU.out
OMP_NUM_THREADS=1
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 4 2 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_2_zins_1CPU.out
./NativeCppBlackScholesSolver solveND 4 3 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_3_zins_1CPU.out
./NativeCppBlackScholesSolver solveND 4 4 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_4_zins_1CPU.out
./NativeCppBlackScholesSolver solveND 4 5 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_5_zins_1CPU.out
./NativeCppBlackScholesSolver solveND 4 6 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_6_zins_1CPU.out
./NativeCppBlackScholesSolver solveND 4 7 BStest4DMC.bound BStest4DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test4D_level_7_zins_1CPU.out

#!/bin/sh

OMP_NUM_THREADS=8
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 2 2 BStest2DMC.bound BStest2DMC_zins.stoch BStest2D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_2_zins_8CPU.out
./NativeCppBlackScholesSolver solveND 2 3 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_3_zins_8CPU.out
./NativeCppBlackScholesSolver solveND 2 4 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_4_zins_8CPU.out
./NativeCppBlackScholesSolver solveND 2 5 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_5_zins_8CPU.out
./NativeCppBlackScholesSolver solveND 2 6 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_6_zins_8CPU.out
./NativeCppBlackScholesSolver solveND 2 7 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_7_zins_8CPU.out
./NativeCppBlackScholesSolver solveND 2 8 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_8_zins_8CPU.out
./NativeCppBlackScholesSolver solveND 2 9 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_9_zins_8CPU.out
./NativeCppBlackScholesSolver solveND 2 10 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_10_zins_8CPU.out
OMP_NUM_THREADS=4
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 2 2 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_2_zins_4CPU.out
./NativeCppBlackScholesSolver solveND 2 3 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_3_zins_4CPU.out
./NativeCppBlackScholesSolver solveND 2 4 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_4_zins_4CPU.out
./NativeCppBlackScholesSolver solveND 2 5 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_5_zins_4CPU.out
./NativeCppBlackScholesSolver solveND 2 6 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_6_zins_4CPU.out
./NativeCppBlackScholesSolver solveND 2 7 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_7_zins_4CPU.out
./NativeCppBlackScholesSolver solveND 2 8 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_8_zins_4CPU.out
./NativeCppBlackScholesSolver solveND 2 9 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_9_zins_4CPU.out
./NativeCppBlackScholesSolver solveND 2 10 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_10_zins_4CPU.out
OMP_NUM_THREADS=2
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 2 2 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_2_zins_2CPU.out
./NativeCppBlackScholesSolver solveND 2 3 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_3_zins_2CPU.out
./NativeCppBlackScholesSolver solveND 2 4 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_4_zins_2CPU.out
./NativeCppBlackScholesSolver solveND 2 5 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_5_zins_2CPU.out
./NativeCppBlackScholesSolver solveND 2 6 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_6_zins_2CPU.out
./NativeCppBlackScholesSolver solveND 2 7 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_7_zins_2CPU.out
./NativeCppBlackScholesSolver solveND 2 8 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_8_zins_2CPU.out
./NativeCppBlackScholesSolver solveND 2 9 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_9_zins_2CPU.out
./NativeCppBlackScholesSolver solveND 2 10 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_10_zins_2CPU.out
OMP_NUM_THREADS=1
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 2 2 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_2_zins_1CPU.out
./NativeCppBlackScholesSolver solveND 2 3 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_3_zins_1CPU.out
./NativeCppBlackScholesSolver solveND 2 4 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_4_zins_1CPU.out
./NativeCppBlackScholesSolver solveND 2 5 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_5_zins_1CPU.out
./NativeCppBlackScholesSolver solveND 2 6 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_6_zins_1CPU.out
./NativeCppBlackScholesSolver solveND 2 7 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_7_zins_1CPU.out
./NativeCppBlackScholesSolver solveND 2 8 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_8_zins_1CPU.out
./NativeCppBlackScholesSolver solveND 2 9 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_9_zins_1CPU.out
./NativeCppBlackScholesSolver solveND 2 10 BStest2DMC.bound BStest2DMC_zins.stoch 1.0 std_euro_call 0.05 1.0 0.1 CrNic 16000 0.00001 > test2D_level_10_zins_1CPU.out
#!/bin/sh

OMP_NUM_THREADS=8
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 2 2 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_2_8CPU.out
./NativeCppBlackScholesSolver solveND 2 3 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_3_8CPU.out
./NativeCppBlackScholesSolver solveND 2 4 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_4_8CPU.out
./NativeCppBlackScholesSolver solveND 2 5 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_5_8CPU.out
./NativeCppBlackScholesSolver solveND 2 6 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_6_8CPU.out
./NativeCppBlackScholesSolver solveND 2 7 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_7_8CPU.out
./NativeCppBlackScholesSolver solveND 2 8 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_8_8CPU.out
./NativeCppBlackScholesSolver solveND 2 9 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_9_8CPU.out
./NativeCppBlackScholesSolver solveND 2 10 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_10_8CPU.out
OMP_NUM_THREADS=4
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 2 2 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_2_4CPU.out
./NativeCppBlackScholesSolver solveND 2 3 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_3_4CPU.out
./NativeCppBlackScholesSolver solveND 2 4 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_4_4CPU.out
./NativeCppBlackScholesSolver solveND 2 5 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_5_4CPU.out
./NativeCppBlackScholesSolver solveND 2 6 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_6_4CPU.out
./NativeCppBlackScholesSolver solveND 2 7 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_7_4CPU.out
./NativeCppBlackScholesSolver solveND 2 8 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_8_4CPU.out
./NativeCppBlackScholesSolver solveND 2 9 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_9_4CPU.out
./NativeCppBlackScholesSolver solveND 2 10 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_10_4CPU.out
OMP_NUM_THREADS=2
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 2 2 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_2_2CPU.out
./NativeCppBlackScholesSolver solveND 2 3 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_3_2CPU.out
./NativeCppBlackScholesSolver solveND 2 4 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_4_2CPU.out
./NativeCppBlackScholesSolver solveND 2 5 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_5_2CPU.out
./NativeCppBlackScholesSolver solveND 2 6 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_6_2CPU.out
./NativeCppBlackScholesSolver solveND 2 7 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_7_2CPU.out
./NativeCppBlackScholesSolver solveND 2 8 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_8_2CPU.out
./NativeCppBlackScholesSolver solveND 2 9 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_9_2CPU.out
./NativeCppBlackScholesSolver solveND 2 10 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_10_2CPU.out
OMP_NUM_THREADS=1
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 2 2 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_2_1CPU.out
./NativeCppBlackScholesSolver solveND 2 3 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_3_1CPU.out
./NativeCppBlackScholesSolver solveND 2 4 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_4_1CPU.out
./NativeCppBlackScholesSolver solveND 2 5 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_5_1CPU.out
./NativeCppBlackScholesSolver solveND 2 6 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_6_1CPU.out
./NativeCppBlackScholesSolver solveND 2 7 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_7_1CPU.out
./NativeCppBlackScholesSolver solveND 2 8 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_8_1CPU.out
./NativeCppBlackScholesSolver solveND 2 9 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_9_1CPU.out
./NativeCppBlackScholesSolver solveND 2 10 BStest2DMC.bound BStest2DMC.stoch BStest2D.strike avgM 0.00 1.0 0.1 CrNic 16000 0.00001 > test2D_level_10_1CPU.out
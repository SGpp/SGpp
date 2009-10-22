#!/bin/sh

OMP_NUM_THREADS=8
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 3 2 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_2_zins_8CPU.out
./NativeCppBlackScholesSolver solveND 3 3 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_3_zins_8CPU.out
./NativeCppBlackScholesSolver solveND 3 4 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_4_zins_8CPU.out
./NativeCppBlackScholesSolver solveND 3 5 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_5_zins_8CPU.out
./NativeCppBlackScholesSolver solveND 3 6 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_6_zins_8CPU.out
./NativeCppBlackScholesSolver solveND 3 7 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_7_zins_8CPU.out
./NativeCppBlackScholesSolver solveND 3 8 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_8_zins_8CPU.out
OMP_NUM_THREADS=4
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 3 2 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_2_zins_4CPU.out
./NativeCppBlackScholesSolver solveND 3 3 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_3_zins_4CPU.out
./NativeCppBlackScholesSolver solveND 3 4 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_4_zins_4CPU.out
./NativeCppBlackScholesSolver solveND 3 5 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_5_zins_4CPU.out
./NativeCppBlackScholesSolver solveND 3 6 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_6_zins_4CPU.out
./NativeCppBlackScholesSolver solveND 3 7 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_7_zins_4CPU.out
./NativeCppBlackScholesSolver solveND 3 8 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_8_zins_4CPU.out
OMP_NUM_THREADS=2
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 3 2 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_2_zins_2CPU.out
./NativeCppBlackScholesSolver solveND 3 3 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_3_zins_2CPU.out
./NativeCppBlackScholesSolver solveND 3 4 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_4_zins_2CPU.out
./NativeCppBlackScholesSolver solveND 3 5 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_5_zins_2CPU.out
./NativeCppBlackScholesSolver solveND 3 6 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_6_zins_2CPU.out
./NativeCppBlackScholesSolver solveND 3 7 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_7_zins_2CPU.out
./NativeCppBlackScholesSolver solveND 3 8 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_8_zins_2CPU.out
OMP_NUM_THREADS=1
export OMP_NUM_THREADS
./NativeCppBlackScholesSolver solveND 3 2 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_2_zins_1CPU.out
./NativeCppBlackScholesSolver solveND 3 3 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_3_zins_1CPU.out
./NativeCppBlackScholesSolver solveND 3 4 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_4_zins_1CPU.out
./NativeCppBlackScholesSolver solveND 3 5 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_5_zins_1CPU.out
./NativeCppBlackScholesSolver solveND 3 6 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_6_zins_1CPU.out
./NativeCppBlackScholesSolver solveND 3 7 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_7_zins_1CPU.out
./NativeCppBlackScholesSolver solveND 3 8 BStest3DMC.bound BStest3DMC_zins.stoch BStest3D.strike avgM 0.05 1.0 0.1 CrNic 16000 0.00001 > test3D_level_8_zins_1CPU.out

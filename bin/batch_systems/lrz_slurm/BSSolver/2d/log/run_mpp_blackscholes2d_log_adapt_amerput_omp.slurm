#!/bin/bash
#SBATCH -o /home/cluster/t1207/lu28mec/SGpp/bin/blackscholes2d_log_adapt_amerput_omp_mpp.%j.out
#SBATCH -D /home/cluster/t1207/lu28mec/SGpp/bin
#SBATCH -J blackscholes2d_log_adapt_amerput_omp_mpp
#SBATCH --partition mpp1_batch
#SBATCH --nodes=1-1
#SBATCH --mail-type=end
#SBATCH --mail-user=heinecke@in.tum.de
#SBATCH --export=NONE
#SBATCH --time=48:00:00

source /etc/profile
cd $HOME/SGpp/bin

module unload ccomp/intel
module load ccomp/intel/12.0

export OMP_NUM_THREADS=16

./BSSolver solveNDadaptSurplus log 2 2 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 2 4 0.000005 coarsen 0.00000001
./BSSolver solveNDadaptSurplus log 2 3 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 3 4 0.000005 coarsen 0.00000001
./BSSolver solveNDadaptSurplus log 2 4 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 4 4 0.000005 coarsen 0.00000001
./BSSolver solveNDadaptSurplus log 2 5 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 5 4 0.000005 coarsen 0.00000001
./BSSolver solveNDadaptSurplus log 2 6 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 6 4 0.000005 coarsen 0.00000001
./BSSolver solveNDadaptSurplus log 2 7 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 7 4 0.000005 coarsen 0.00000001
./BSSolver solveNDadaptSurplus log 2 8 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 8 4 0.000005 coarsen 0.00000001
./BSSolver solveNDadaptSurplus log 2 9 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 9 4 0.000005 coarsen 0.00000001
./BSSolver solveNDadaptSurplus log 2 10 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 10 4 0.000005 coarsen 0.00000001
./BSSolver solveNDadaptSurplus log 2 11 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 11 4 0.000005 coarsen 0.00000001

./BSSolver solveNDadaptSurplusSubDomain log 2 2 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 2 4 0.000005 coarsen 0.00000001
./BSSolver solveNDadaptSurplusSubDomain log 2 3 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 3 4 0.000005 coarsen 0.00000001
./BSSolver solveNDadaptSurplusSubDomain log 2 4 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 4 4 0.000005 coarsen 0.00000001
./BSSolver solveNDadaptSurplusSubDomain log 2 5 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 5 4 0.000005 coarsen 0.00000001
./BSSolver solveNDadaptSurplusSubDomain log 2 6 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 6 4 0.000005 coarsen 0.00000001
./BSSolver solveNDadaptSurplusSubDomain log 2 7 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 7 4 0.000005 coarsen 0.00000001
./BSSolver solveNDadaptSurplusSubDomain log 2 8 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 8 4 0.000005 coarsen 0.00000001
./BSSolver solveNDadaptSurplusSubDomain log 2 9 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 9 4 0.000005 coarsen 0.00000001
./BSSolver solveNDadaptSurplusSubDomain log 2 10 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 10 4 0.000005 coarsen 0.00000001
./BSSolver solveNDadaptSurplusSubDomain log 2 11 BStest2d.log.bound BStest2d.stoch 1.0 std_amer_put 0.05 1.0 0.025 CrNic 8000 0.000001 maxLevel 11 4 0.000005 coarsen 0.00000001
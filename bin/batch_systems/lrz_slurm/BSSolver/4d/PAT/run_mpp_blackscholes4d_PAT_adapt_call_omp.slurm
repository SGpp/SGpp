#!/bin/bash
#SBATCH -o /home/cluster/t1207/lu28mec/SGpp/bin/blackscholes4d_PAT_adapt_call_omp_mpp.%j.out
#SBATCH -D /home/cluster/t1207/lu28mec/SGpp_DM/bin
#SBATCH -J blackscholes4d_PAT_adaptdomain_call_omp_mpp
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

./BSSolver solveNDadaptSurplus PAT 4 2 BStest4d.log.bound BStest4d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 160 0.000001 maxLevel 2 5 1e-6 none 1e-7
./BSSolver solveNDadaptSurplus PAT 4 3 BStest4d.log.bound BStest4d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 160 0.000001 maxLevel 3 5 1e-6 none 1e-7
./BSSolver solveNDadaptSurplus PAT 4 4 BStest4d.log.bound BStest4d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 160 0.000001 maxLevel 4 5 1e-6 none 1e-7
./BSSolver solveNDadaptSurplus PAT 4 5 BStest4d.log.bound BStest4d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 160 0.000001 maxLevel 5 5 1e-6 none 1e-7
./BSSolver solveNDadaptSurplus PAT 4 6 BStest4d.log.bound BStest4d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 160 0.000001 maxLevel 6 5 1e-6 none 1e-7

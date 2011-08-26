#!/bin/bash
#SBATCH -o /home/cluster/t1207/lu28mec/SGpp/bin/blackscholes3d_PAT_adapt_put_omp_mpp.%j.out
#SBATCH -D /home/cluster/t1207/lu28mec/SGpp/bin
#SBATCH -J blackscholes3d_PAT_adapt_put_omp_mpp
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

./BSSolver solveNDadaptSurplus PAT 3 2 BStest3d.log.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 8000 0.000001 maxLevel 2 4 0.000001 coarsen 0.0000001
./BSSolver solveNDadaptSurplus PAT 3 3 BStest3d.log.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 8000 0.000001 maxLevel 3 4 0.000001 coarsen 0.0000001
./BSSolver solveNDadaptSurplus PAT 3 4 BStest3d.log.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 8000 0.000001 maxLevel 4 4 0.000001 coarsen 0.0000001
./BSSolver solveNDadaptSurplus PAT 3 5 BStest3d.log.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 8000 0.000001 maxLevel 5 4 0.000001 coarsen 0.0000001
./BSSolver solveNDadaptSurplus PAT 3 6 BStest3d.log.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 8000 0.000001 maxLevel 6 4 0.000001 coarsen 0.0000001
./BSSolver solveNDadaptSurplus PAT 3 7 BStest3d.log.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 8000 0.000001 maxLevel 7 4 0.000001 coarsen 0.0000001
./BSSolver solveNDadaptSurplus PAT 3 8 BStest3d.log.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 8000 0.000001 maxLevel 8 4 0.000001 coarsen 0.0000001
./BSSolver solveNDadaptSurplus PAT 3 9 BStest3d.log.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 8000 0.000001 maxLevel 9 4 0.000001 coarsen 0.0000001

./BSSolver solveNDadaptSurplusSubDomain PAT 3 2 BStest3d.log.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 20000 0.000001 maxLevel 2 4 0.000001 coarsen 0.0000001
./BSSolver solveNDadaptSurplusSubDomain PAT 3 3 BStest3d.log.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 20000 0.000001 maxLevel 3 4 0.000001 coarsen 0.0000001
./BSSolver solveNDadaptSurplusSubDomain PAT 3 4 BStest3d.log.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 20000 0.000001 maxLevel 4 4 0.000001 coarsen 0.0000001
./BSSolver solveNDadaptSurplusSubDomain PAT 3 5 BStest3d.log.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 20000 0.000001 maxLevel 5 4 0.000001 coarsen 0.0000001
./BSSolver solveNDadaptSurplusSubDomain PAT 3 6 BStest3d.log.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 20000 0.000001 maxLevel 6 4 0.000001 coarsen 0.0000001
./BSSolver solveNDadaptSurplusSubDomain PAT 3 7 BStest3d.log.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 20000 0.000001 maxLevel 7 4 0.000001 coarsen 0.0000001
./BSSolver solveNDadaptSurplusSubDomain PAT 3 8 BStest3d.log.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 20000 0.000001 maxLevel 8 4 0.000001 coarsen 0.0000001
./BSSolver solveNDadaptSurplusSubDomain PAT 3 9 BStest3d.log.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 20000 0.000001 maxLevel 9 4 0.000001 coarsen 0.0000001
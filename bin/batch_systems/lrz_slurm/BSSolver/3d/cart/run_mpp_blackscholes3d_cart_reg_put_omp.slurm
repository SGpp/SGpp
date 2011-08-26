#!/bin/bash
#SBATCH -o /home/cluster/t1207/lu28mec/SGpp/bin/blackscholes3d_cart_reg_put_omp_mpp.%j.out
#SBATCH -D /home/cluster/t1207/lu28mec/SGpp/bin
#SBATCH -J blackscholes3d_cart_reg_put_omp_mpp
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

./BSSolver solveNDanalyze cart 3 3 13 BStest3d.cart.bound BStest3d.stoch 1.0 std_euro_put 0.05 1.0 0.05 CrNic 8000 0.000001 BStest3d.analyze

#!/bin/bash
#SBATCH -o /home/cluster/t1207/lu28mec/SGpp/bin/blackscholes2d_cart_reg_call_omp_mpp.%j.out
#SBATCH -D /home/cluster/t1207/lu28mec/SGpp/bin
#SBATCH -J blackscholes2d_cart_reg_call_omp_mpp
#SBATCH --partition mpp1_batch
#SBATCH --nodes=1-1
#SBATCH --mail-type=end
#SBATCH --mail-user=heinecke@in.tum.de
#SBATCH --export=NONE
#SBATCH --time=48:00:00

source /etc/profile
cd $HOME/SGpp/bin

module unload ccomp/intel
module load ccomp/intel/12.1

export OMP_NUM_THREADS=16

./BSSolver solveNDanalyze cart 2 3 14 BStest2d.cart.bound BStest2d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 8000 0.000001 BStest2d.analyze

#!/bin/bash
#SBATCH -o /home/cluster/t1207/lu28mec/SGpp/bin/blackscholes6d_PAT_adaptdomain_amerput_mpi_7r_16thr_mpp.%j.out
#SBATCH -D /home/cluster/t1207/lu28mec/SGpp_DM/bin
#SBATCH -J blackscholes6d_PAT_adaptdomain_amerput_mpi_7r_16thr_mpp
#SBATCH --partition mpp1_batch
#SBATCH --ntasks=7
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=end
#SBATCH --mail-user=heinecke@in.tum.de
#SBATCH --export=NONE
#SBATCH --time=48:00:00

source /etc/profile
cd $HOME/SGpp/bin

module unload ccomp/intel
module load ccomp/intel/12.0

export OMP_NUM_THREADS=16
export I_MPI_PIN_MODE=mpd
export I_MPI_DEBUG=4
export PSI_EXPORTS=$PSI_EXPORTS,I_MPI_PIN_MODE,I_MPI_DEBUG

srun_ps -v -t 16 ./BSSolverMPI solveNDadaptSurplusSubDomain PAT 6 2 BStest6d.log.bound BStest6d.stoch 1.0 std_amer_put 0.05 1.0 0.0125 CrNic 180 0.000001 maxLevel 2 8 0.000001 none 0.0000001
srun_ps -v -t 16 ./BSSolverMPI solveNDadaptSurplusSubDomain PAT 6 3 BStest6d.log.bound BStest6d.stoch 1.0 std_amer_put 0.05 1.0 0.0125 CrNic 180 0.000001 maxLevel 3 8 0.000001 none 0.0000001
srun_ps -v -t 16 ./BSSolverMPI solveNDadaptSurplusSubDomain PAT 6 4 BStest6d.log.bound BStest6d.stoch 1.0 std_amer_put 0.05 1.0 0.0125 CrNic 180 0.000001 maxLevel 4 8 0.000001 none 0.0000001
srun_ps -v -t 16 ./BSSolverMPI solveNDadaptSurplusSubDomain PAT 6 5 BStest6d.log.bound BStest6d.stoch 1.0 std_amer_put 0.05 1.0 0.0125 CrNic 180 0.000001 maxLevel 5 8 0.000001 none 0.0000001
srun_ps -v -t 16 ./BSSolverMPI solveNDadaptSurplusSubDomain PAT 6 6 BStest6d.log.bound BStest6d.stoch 1.0 std_amer_put 0.05 1.0 0.0125 CrNic 180 0.000001 maxLevel 6 8 0.000001 none 0.0000001

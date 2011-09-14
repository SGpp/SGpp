#!/bin/bash
#SBATCH -o /home/cluster/t1207/lu28mec/SGpp/bin/blackscholes5d_log_adaptdomain_call_mpi_6r_08thr_mpp.%j.out
#SBATCH -D /home/cluster/t1207/lu28mec/SGpp_DM/bin
#SBATCH -J blackscholes5d_log_adaptdomain_call_mpi_6r_08thr_mpp
#SBATCH --partition mpp1_batch
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=end
#SBATCH --mail-user=heinecke@in.tum.de
#SBATCH --export=NONE
#SBATCH --time=48:00:00

source /etc/profile
cd $HOME/SGpp/bin

module unload ccomp/intel
module load ccomp/intel/12.1

export OMP_NUM_THREADS=8
export I_MPI_PIN_MODE=mpd
export I_MPI_DEBUG=4
export PSI_EXPORTS=$PSI_EXPORTS,I_MPI_PIN_MODE,I_MPI_DEBUG

srun_ps -v -t 8 ./BSSolverMPI solveNDadaptSurplusSubDomain log 5 2 BStest5d.log.bound BStest5d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 180 0.000001 maxLevel 3 7 0.000001 coarsen 0.0000001
srun_ps -v -t 8 ./BSSolverMPI solveNDadaptSurplusSubDomain log 5 3 BStest5d.log.bound BStest5d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 180 0.000001 maxLevel 4 7 0.000001 coarsen 0.0000001
srun_ps -v -t 8 ./BSSolverMPI solveNDadaptSurplusSubDomain log 5 4 BStest5d.log.bound BStest5d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 180 0.000001 maxLevel 5 7 0.000001 coarsen 0.0000001
srun_ps -v -t 8 ./BSSolverMPI solveNDadaptSurplusSubDomain log 5 5 BStest5d.log.bound BStest5d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 180 0.000001 maxLevel 6 7 0.000001 coarsen 0.0000001
srun_ps -v -t 8 ./BSSolverMPI solveNDadaptSurplusSubDomain log 5 6 BStest5d.log.bound BStest5d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 180 0.000001 maxLevel 7 7 0.000001 coarsen 0.0000001
srun_ps -v -t 8 ./BSSolverMPI solveNDadaptSurplusSubDomain log 5 7 BStest5d.log.bound BStest5d.stoch 1.0 std_euro_call 0.05 1.0 0.05 CrNic 180 0.000001 maxLevel 8 7 0.000001 coarsen 0.0000001

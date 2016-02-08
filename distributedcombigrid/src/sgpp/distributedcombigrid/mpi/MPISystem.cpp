/*
 * MPISetup.cpp
 *
 *  Created on: Jan 23, 2013
 *      Author: mh
 *
 *  Partially copied from the pe Physics Engine class MPISystem
 */

#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"

#include <stdexcept>

namespace combigrid {

/*!\brief Constructor for the MPISystem class.
 //
 // \exception std::runtime_error MPI was not initialized
 //
 // Constructor for the MPI System. The default global communicator and local communicator is MPI_COMM_WORLD.
 // The total number of MPI processes and the rank of the MPI process in is determined from
 // the communicator. The default root process has rank zero.
 */
MPISystem::MPISystem() {
  // check if MPI was initialized (e.g. by MPI_Init or similar)
  int mpiInitialized(0);
  MPI_Initialized(&mpiInitialized);

  if (!mpiInitialized)
    throw std::runtime_error("MPI was not initialized");

  // Setting the communicator
  globalComm_ = MPI_COMM_WORLD;

  // Reinitializing the total number of processes and the rank of this process
  MPI_Comm_size(globalComm_,
                &globalSize_); // Setting the total number of MPI processes
  MPI_Comm_rank(globalComm_, &globalRank_);  // Setting the rank of this process

  // Setting the root process to rank 0
  globalRoot_ = 0;

  // Setting the local communicator to the same settings
  localComm_ = globalComm_;
  localSize_ = globalSize_;
  localRank_ = globalRank_;
  localRoot_ = globalRoot_;

  active_.resize(globalSize_, true);
}

/*!\brief Destructor for the MPISystem class.
 */
MPISystem::~MPISystem() {
}

/*!\brief Output of the current state of the MPI system.
 */
void MPISystem::print(std::ostream& os) const {
  os << " Number of MPI processes in global communicator = " << globalSize_
     << "\n" << " My rank in global communicator = " << globalRank_ << "\n"
     << " \n" << " Number of MPI processes in local communicator = "
     << localSize_ << "\n" << " My rank in local communicator = " << localRank_
     << "\n";
}

/*!\brief Sets the root process in the global communicator
 //
 // \param rootProcess The new root process in the global communicator.
 // \return void
 // \exception std::invalid_argument Invalid root process.
 //
 // This function sets the root processes of the MPI system. In case \a rootProcess
 // is larger or equal than the total number of processes, a \a std::invalid_argument exception
 // is thrown.
 */
void MPISystem::setGlobalRoot(RankType rootProcess) {
  if (rootProcess < 0 || rootProcess >= globalSize_)
    throw std::invalid_argument("Invalid root process");

  globalRoot_ = rootProcess;
}

/*!\brief Setting the global MPI communicator of the MPI system.
 //
 // \param communicator The new global MPI communicator.
 // \return void
 //
 // This function sets the global MPI communicator of the MPI parallel simulation. Additionally,
 // it reevaluates the total number of processes and the rank of this processes for the
 // new communicator. The root process in the global communicator is reset to zero.
 */
void MPISystem::setGlobalComm(CommunicatorType communicator) {
  // Setting the communicator
  globalComm_ = communicator;

  // Reinitializing the total number of processes and the rank of this process
  MPI_Comm_size(globalComm_,
                &globalSize_); // Estimating the total number of MPI processes
  MPI_Comm_rank(globalComm_, &globalRank_); // Estimating the rank of this process

  // Reset root process to zero
  globalRoot_ = 0;
}

/*!\brief Sets the root process in the local communicator
 //
 // \param rootProcess The new root process in the local communicator.
 // \return void
 // \exception std::invalid_argument Invalid root process.
 //
 // This function sets the root processes of the MPI system. In case \a rootProcess
 // is larger or equal than the total number of processes, a \a std::invalid_argument exception
 // is thrown.
 */
void MPISystem::setLocalRoot(RankType rootProcess) {
  if (rootProcess < 0 || rootProcess >= localSize_)
    throw std::invalid_argument("Invalid root process");

  localRoot_ = rootProcess;
}

/*!\brief Setting the local MPI communicator of the MPI system.
 //
 // \param communicator The new local MPI communicator.
 // \return void
 //
 // This function sets the global MPI communicator of the MPI parallel simulation. Additionally,
 // it reevaluates the total number of processes and the rank of this processes for the
 // new communicator. The root process in the global communicator is reset to zero.
 */
void MPISystem::setLocalComm(CommunicatorType communicator) {
  // Setting the communicator
  localComm_ = communicator;

  // Reinitializing the total number of processes and the rank of this process
  MPI_Comm_size(localComm_,
                &localSize_); // Estimating the total number of MPI processes
  MPI_Comm_rank(localComm_, &localRank_); // Estimating the rank of this process

  // Reset root process to zero
  localRoot_ = 0;
}

/*!\brief Global output operator for the MPI system.
 // \ingroup mpi
 //
 // \param os Reference to the output stream.
 // \param ms Reference to a constant MPI system object.
 // \return Reference to the output stream.
 */
std::ostream& operator<<(std::ostream& os, const MPISystem& ms) {
  os << "--" << "MPI SYSTEM PARAMETERS"
     << "---------------------------------------------------------\n";
  ms.print(os);
  os
      << "--------------------------------------------------------------------------------\n"
      << std::endl;
  return os;
}

/*!\brief Global output operator for MPI system handles.
 // \ingroup mpi
 //
 // \param os Reference to the output stream.
 // \param ms MPI system handle.
 // \return Reference to the output stream.
 */
std::ostream& operator<<(std::ostream& os, const MPISystemID& ms) {
  os << "--" << "MPI SYSTEM PARAMETERS"
     << "---------------------------------------------------------\n";
  ms->print(os);
  os
      << "--------------------------------------------------------------------------------\n"
      << std::endl;
  return os;
}

/*!\brief Global output operator for constant MPI system handles.
 // \ingroup mpi
 //
 // \param os Reference to the output stream.
 // \param ms Constant MPI system handle.
 // \return Reference to the output stream.
 */
std::ostream& operator<<(std::ostream& os, const ConstMPISystemID& ms) {
  os << "--" << "MPI SYSTEM PARAMETERS"
     << "---------------------------------------------------------\n";
  ms->print(os);
  os
      << "--------------------------------------------------------------------------------\n"
      << std::endl;
  return os;
}

} // namespace combigrid

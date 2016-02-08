#ifndef MPISYSTEM_HPP
#define MPISYSTEM_HPP

#include <mpi.h>
#include <ostream>
#include <vector>
#include <assert.h>

#include "sgpp/distributedcombigrid/mpi/MPISystemID.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

namespace combigrid {

/** MPI Communication System
 *
 * the class is partially take from the pe physics engine
 *
 * The MPISystem class represents the MPI communication system.
 */
class MPISystem {
 private:
  /** Constructor
   *
   */
  explicit MPISystem();

 public:
  /** Destructor
   *
   */
  virtual ~MPISystem();

  inline CommunicatorType getGlobalComm() const;
  inline int getGlobalSize() const;
  inline RankType getGlobalID() const;
  inline RankType getManagerID() const;

  inline CommunicatorType getLocalComm() const;
  inline int getLocalSize() const;
  inline RankType getLocalID() const;
  inline RankType getLocalRootID() const;

  void setGlobalRoot(RankType rootID);
  void setGlobalComm(CommunicatorType comm);

  void setLocalRoot(RankType localRootID);
  void setLocalComm(CommunicatorType comm);

  /** Print MPI System information
   *
   * @param os output stream where the mpi system information shall be printed
   */
  void print(std::ostream& os) const;
  //@}

  inline bool checkActive(RankType rank) const;
 private:
  friend MPISystemID theMPISystem();

  CommunicatorType globalComm_;   //!< The global MPI communicator
  int globalSize_;   //!< Number of MPI processes in the global communicator
  RankType globalRank_; //!< Rank of the current MPI process in the global communicator
  RankType globalRoot_;   //!< Root MPI process in the global communicator

  CommunicatorType localComm_;   //!< The local MPI communicator
  int localSize_;   //!< Number of MPI processes in the local communicator
  RankType localRank_; //!< Rank of the current MPI process in the local communicator
  RankType localRoot_;   //!< Root MPI process in the local communicator

  std::vector<bool> active_;
};

/*!\name MPI communication system setup functions */
//@{
inline MPISystemID theMPISystem();
//@}

/** Returns a handle to the MPI communication system.
 *
 * This function returns a handle to the MPI communication system. This handle
 * can be used to configure the communication system or to acquire
 * the current settings. The function expects that MPI has already been properly
 * initialized (e.g. via the MPI_Init() or any similar function). In case MPI
 * was not initialized, a \a std::runtime_error exception is thrown.
 *
 * \return Handle to the MPI communication system.
 * \exception std::runtime_error MPI system is not initialized.
 */
inline MPISystemID theMPISystem() {
  static MPISystemID system(new MPISystem());
  return system;
}

/** Returns the current global MPI communicator
 *
 */
inline CommunicatorType MPISystem::getGlobalComm() const {
  return globalComm_;
}

/** Returns the total number of MPI processes in the global communicator.
 *
 */
inline int MPISystem::getGlobalSize() const {
  return globalSize_;
}

/** Returns the rank of the current MPI process in the global communicator.
 *
 */
inline RankType MPISystem::getGlobalID() const {
  return globalRank_;
}

/** Returns the root process in the global communicator
 *
 */
inline RankType MPISystem::getManagerID() const {
  return globalRoot_;
}

/** Returns the current local MPI communicator
 *
 */
inline CommunicatorType MPISystem::getLocalComm() const {
  return localComm_;
}

/** Returns the total number of MPI processes in the local communicator.
 *
 */
inline int MPISystem::getLocalSize() const {
  return localSize_;
}

/** Returns the rank of the current MPI process in the local communicator.
 *
 */
inline RankType MPISystem::getLocalID() const {
  return localRank_;
}

/** Returns the rank of the current MPI process in the local communicator.
 *
 */
inline RankType MPISystem::getLocalRootID() const {
  return localRoot_;
}

inline bool MPISystem::checkActive(int rank) const {
  assert(rank >= 0);
  assert(rank < globalSize_);
  return active_[rank];
}

// operators
std::ostream& operator<<(std::ostream& os, const MPISystem& ms);
std::ostream& operator<<(std::ostream& os, const MPISystemID& ms);
std::ostream& operator<<(std::ostream& os, const ConstMPISystemID& ms);

} //namespace combigrid

#endif // MPISYSTEM_HPP

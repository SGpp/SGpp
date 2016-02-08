/*
 * ProcessGroupManager.hpp
 *
 *  Created on: Jul 17, 2014
 *      Author: heenemo
 */

#ifndef PROCESSGROUPMANAGER_HPP_
#define PROCESSGROUPMANAGER_HPP_

#include <vector>

#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/combicom/CombiCom.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupSignals.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"

namespace combigrid {

class ProcessGroupManager {

 public:
  ProcessGroupManager(RankType managerID, RankType pgroupRootID,
                      CommunicatorType globalComm);

  virtual
  ~ProcessGroupManager();

  bool
  runfirst(Task* t);

  bool
  runnext();

  bool
  exit();

  inline StatusType
  getStatus();

  inline complex
  eval(const std::vector<real>& coords);

  bool
  sync();

  inline const TaskContainer&
  getTaskContainer() const;

  bool
  combine();

  template<typename FG_ELEMENT>
  bool
  combineFG(FullGrid<FG_ELEMENT>& fg);

  template<typename FG_ELEMENT>
  bool
  gridEval(FullGrid<FG_ELEMENT>& fg);

  bool
  gridGather(LevelVector& leval);

  bool
  updateCombiParameters(CombiParameters& params);

  bool addTask( Task* );

  bool recompute( Task* );

 private:

  RankType managerID_; // rank of manager

  RankType pgroupRootID_; // rank of the root process of the pgroup

  CommunicatorType gcomm_; // global communicator which includes at least
  // manager and pgroup roots

  TaskContainer tasks_;

  StatusType status_;

  MPI_Request statusRequest_;

  //TaskIDContainer EVtasks_;
};

typedef std::vector<ProcessGroupManager> ProcessGroupManagerContainer;

inline StatusType ProcessGroupManager::getStatus() {
  if (status_ == PROCESS_GROUP_WAIT)
    return PROCESS_GROUP_WAIT;

  /* todo: do we really need the exit status?
   if( status_ == APP_EXIT )
   return APP_EXIT;
   */

  // todo: look into details of the mpi standard whether we need the function
  // mpi_test, like mpi_wait, completes a non-blocking communication
  // it is probably not necessary, but might give mpi the chance to free some
  // internally assigned buffers
  // it also might avoid getting stuck in infinite loops when repeatedly
  // checking status without a system call
  int flag;
  MPI_Test(&statusRequest_, &flag, MPI_STATUS_IGNORE);

  return status_;
}

inline complex ProcessGroupManager::eval(const std::vector<real>& x) {
  //todo: implement
  return complex(0.0, 0.0);
}

inline const TaskContainer&
ProcessGroupManager::getTaskContainer() const {
  return tasks_;
}

template<typename FG_ELEMENT>
bool ProcessGroupManager::gridEval(FullGrid<FG_ELEMENT>& fg) {
  // can only send sync signal when in wait state, so check first
  assert(status_ == PROCESS_GROUP_WAIT);

  // send signal
  SignalType signal = GRID_EVAL;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, gcomm_);

  // send levelvector
  std::vector<int> tmp(fg.getLevels().begin(), fg.getLevels().end());
  MPI_Send(&tmp[0], static_cast<int>(tmp.size()), MPI_INT, pgroupRootID_, 0,
           gcomm_);

  return true;
}

template<typename FG_ELEMENT>
bool ProcessGroupManager::combineFG(FullGrid<FG_ELEMENT>& fg) {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  SignalType signal = COMBINE_FG;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, gcomm_);

  // send levelvector
  std::vector<int>& tmp = fg.getLevels();
  MPI_Send(&tmp[0], static_cast<int>(tmp.size()), MPI_INT, pgroupRootID_, 0,
           gcomm_);

  return true;
}

inline
bool ProcessGroupManager::gridGather(LevelVector& leval) {
  // can only send sync signal when in wait state, so check first
  assert(status_ == PROCESS_GROUP_WAIT);

  // send signal
  SignalType signal = GRID_GATHER;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, gcomm_);

  // send levelvector
  std::vector<int> tmp(leval.begin(), leval.end());
  MPI_Send(&tmp[0], static_cast<int>(tmp.size()), MPI_INT, pgroupRootID_, 0,
           gcomm_);

  return true;
}
} /* namespace combigrid */

#endif /* PROCESSGROUPMANAGER_HPP_ */

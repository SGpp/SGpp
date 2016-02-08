/*
 * ProcessGroupManager.cpp
 *
 *  Created on: Jul 17, 2014
 *      Author: heenemo
 */

#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/mpi/MPIUtils.hpp"

namespace combigrid {
ProcessGroupManager::ProcessGroupManager(RankType managerID,
    RankType pgroupRootID, CommunicatorType globalComm) :
    managerID_(managerID), pgroupRootID_(pgroupRootID), gcomm_(globalComm), status_(
        PROCESS_GROUP_WAIT), statusRequest_(MPI_Request()) {
}

ProcessGroupManager::~ProcessGroupManager() {
  // TODO Auto-generated destructor stub
}

bool ProcessGroupManager::runfirst(Task* t) {
  // first check status
  // tying to add a task to a busy group is an invalid operation
  // and should be avoided
  if (status_ != PROCESS_GROUP_WAIT)
    return false;

  // add task to list of tasks managed by this pgroup
  tasks_.push_back(t);

  // send runfirst_signal to pgroup
  SignalType signal = RUN_FIRST;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, gcomm_);

  // send task
  Task::send(&t, pgroupRootID_, gcomm_);

  // set status
  status_ = PROCESS_GROUP_BUSY;

  // start non-blocking call to receive status
  MPI_Irecv(&status_, 1, MPI_INT, pgroupRootID_, statusTag, gcomm_,
      &statusRequest_);

  // only return true if task successfully send to pgroup
  return true;
}

bool ProcessGroupManager::runnext() {
  // first check status
  // trying to send a command to a busy group is an invalid operation
  // and should be avoided
  assert(status_ == PROCESS_GROUP_WAIT);

  if (tasks_.size() == 0)
    return false;

  // send runnext signal
  SignalType signal = RUN_NEXT;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, gcomm_);

  status_ = PROCESS_GROUP_BUSY;

  // start non-blocking call to receive status
  MPI_Irecv(&status_, 1, MPI_INT, pgroupRootID_, statusTag, gcomm_,
      &statusRequest_);

  return true;
}

bool ProcessGroupManager::exit() {
  // can only send exit signal when in wait state
  if (status_ != PROCESS_GROUP_WAIT)
    return false;

  // send exit signal
  SignalType signal = EXIT;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, gcomm_);

  return true;
}

bool ProcessGroupManager::sync() {
  // can only send sync signal when in wait state
  if (status_ != PROCESS_GROUP_WAIT)
    return false;

  // send sync signal
  SignalType signal = SYNC_TASKS;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, gcomm_);

  for (size_t i = 0; i < tasks_.size(); ++i) {
    Task* t = tasks_[i];

    Task::receive(&t, pgroupRootID_, gcomm_);
  }

  return true;
}

bool ProcessGroupManager::combine() {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  SignalType signal = COMBINE;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, gcomm_);

  return true;
}

bool ProcessGroupManager::updateCombiParameters(CombiParameters& params) {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  SignalType signal = UPDATE_COMBI_PARAMETERS;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, gcomm_);

  // send combiparameters
  MPIUtils::sendClass(&params, pgroupRootID_, gcomm_);

  return true;
}

bool ProcessGroupManager::addTask( Task* t ){
  // first check status
  // tying to add a task to a busy group is an invalid operation
  // and should be avoided
  if (status_ != PROCESS_GROUP_WAIT)
    return false;

  // add task to list of tasks managed by this pgroup
  tasks_.push_back(t);

  // send add task signal to pgroup
  SignalType signal = ADD_TASK;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, gcomm_);

  // send task
  Task::send(&t, pgroupRootID_, gcomm_);

  // set status
  status_ = PROCESS_GROUP_BUSY;

  // start non-blocking call to receive status
  MPI_Irecv(&status_, 1, MPI_INT, pgroupRootID_, statusTag, gcomm_,
      &statusRequest_);

  // only return true if task successfully send to pgroup
  return true;
}


bool ProcessGroupManager::recompute( Task* t ){
  // first check status
  // tying to add a task to a busy group is an invalid operation
  // and should be avoided
  if (status_ != PROCESS_GROUP_WAIT)
    return false;

  // add task to list of tasks managed by this pgroup
  tasks_.push_back(t);

  // send add task signal to pgroup
  SignalType signal = RECOMPUTE;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, gcomm_);

  // send task
  Task::send(&t, pgroupRootID_, gcomm_);

  // set status
  status_ = PROCESS_GROUP_BUSY;

  // start non-blocking call to receive status
  MPI_Irecv(&status_, 1, MPI_INT, pgroupRootID_, statusTag, gcomm_,
      &statusRequest_);

  // only return true if task successfully send to pgroup
  return true;
}

} /* namespace combigrid */

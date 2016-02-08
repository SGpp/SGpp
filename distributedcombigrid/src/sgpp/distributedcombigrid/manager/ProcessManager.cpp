/*
 * ProcessManager.cpp
 *
 *  Created on: Oct 8, 2013
 *      Author: heenemo
 */

#include <algorithm>
#include <iostream>
#include "sgpp/distributedcombigrid/manager/ProcessManager.hpp"
#include "sgpp/distributedcombigrid/combicom/CombiCom.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

namespace combigrid {

bool compareInstances(const Task* instance1, const Task* instance2) {
  return (instance1->estimateRuntime() > instance2->estimateRuntime());
}

ProcessManager::ProcessManager(ProcessGroupManagerContainer& pgroups,
                               TaskContainer& tasks, CombiParameters& params, RankType managerID,
                               CommunicatorType globalComm) :
  pgroups_(pgroups), tasks_(tasks), params_(params), managerID_(managerID),
  gcomm_(
    globalComm) {
}

ProcessManager::~ProcessManager() {
}

void ProcessManager::runfirst() {
  // sort instances in decreasing order
  std::sort(tasks_.begin(), tasks_.end(), compareInstances);

  for (size_t i = 0; i < tasks_.size(); ++i) {
    // wait for available process group
    ProcessGroupManager* g = wait();

    // assign instance to group
    g->runfirst(tasks_[i]);
  }

  std::cout << "ALL TASKS ASSIGNED" << std::endl;

  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i].getStatus() == PROCESS_GROUP_WAIT)
        ++numWaiting;
    }

  }

  std::cout << "RUN_FIRST: all pgroups finished" << std::endl;
}

void ProcessManager::runnext() {
  for (size_t i = 0; i < pgroups_.size(); ++i) {
    pgroups_[i].runnext();
  }

  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i].getStatus() == PROCESS_GROUP_WAIT)
        ++numWaiting;
    }
  }

  std::cout << "RUN_NEXT: all pgroups finished" << std::endl;
}

void ProcessManager::exit() {
  // wait until all process groups are in wait state
  // after sending the exit signal checking the status might not be possible
  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i].getStatus() == PROCESS_GROUP_WAIT)
        ++numWaiting;
    }
  }

  // send exit signal to each group
  for (size_t i = 0; i < pgroups_.size(); ++i) {
    assert(pgroups_[i].exit());
  }
}

void ProcessManager::sync() {
  // container to mark synced pgroups
  std::vector<bool> isGroupSynced(pgroups_.size(), false);

  size_t numSynced(0);

  while (numSynced != pgroups_.size()) {
    numSynced = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      // check only those which are not yet synced
      if (!isGroupSynced[i]) {
        // check status first -> infinity loop avoidance stuff
        if (pgroups_[i].getStatus() == PROCESS_GROUP_WAIT) {
          if (pgroups_[i].sync()) {
            isGroupSynced[i] = true;
            numSynced++;
          }
        }
      }
    }
  }

  std::cout << "All Tasks synced!" << std::endl;
}

void ProcessManager::updateCombiParameters() {
  // container to mark synced pgroups
  std::vector<bool> isGroupSynced(pgroups_.size(), false);

  size_t numSynced(0);

  while (numSynced != pgroups_.size()) {
    numSynced = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      // check only those which are not yet synced
      if (!isGroupSynced[i]) {
        // check status first -> infinity loop avoidance stuff
        if (pgroups_[i].getStatus() == PROCESS_GROUP_WAIT) {
          if (pgroups_[i].updateCombiParameters(params_)) {
            isGroupSynced[i] = true;
            numSynced++;
          }
        }
      }
    }
  }

  std::cout << "All combi parameters updated!" << std::endl;
}


void ProcessManager::redistribute( std::vector<int>& taskID ) {
  for (size_t i = 0; i < taskID.size(); ++i) {
    // find id in list of tasks
    Task* t = NULL;

    for ( Task* tmp : tasks_ ) {
      if ( tmp->getID() == taskID[i] ) {
        t = tmp;
        break;
      }
    }

    assert( t != NULL );

    // wait for available process group
    ProcessGroupManager* g = wait();

    // assign instance to group
    g->addTask( t );
  }

  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i].getStatus() == PROCESS_GROUP_WAIT)
        ++numWaiting;
    }

  }

  std::cout << "Redistribute finished" << std::endl;
}


void ProcessManager::recompute( std::vector<int>& taskID ) {
  for (size_t i = 0; i < taskID.size(); ++i) {
    // find id in list of tasks
    Task* t = NULL;

    for ( Task* tmp : tasks_ ) {
      if ( tmp->getID() == taskID[i] ) {
        t = tmp;
        break;
      }
    }

    assert( t != NULL );

    // wait for available process group
    ProcessGroupManager* g = wait();

    // assign instance to group
    g->recompute( t );
  }

  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i].getStatus() == PROCESS_GROUP_WAIT)
        ++numWaiting;
    }

  }

  std::cout << "Recompute finished" << std::endl;
}

} /* namespace combigrid */


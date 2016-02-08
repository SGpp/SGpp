/*
 * ProcessGroupWorker.hpp
 *
 *  Created on: Jun 24, 2014
 *      Author: heenemo
 */

#ifndef PROCESSGROUPWORKER_HPP_
#define PROCESSGROUPWORKER_HPP_

#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupSignals.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"

namespace combigrid {

class ProcessGroupWorker {
public:
  ProcessGroupWorker(MPI_Comm gcomm, MPI_Comm lcomm, RankType grank,
      RankType lrank, RankType manager, RankType lroot = 0);

  virtual ~ProcessGroupWorker();

  // wait for command from manager
  SignalType wait();

  // send ready signal to manager
  void ready();

  // todo: maybe only needed for gene?
  inline Task* getCurrentTask();

  void combine();

  void combineUniform();

  void combineFG();

  void gridEval();

  void updateCombiParameters();

  /*

   inline CommunicatorType getGlobalComm() const;

   inline CommunicatorType getLocalComm() const;

   inline RankType getGlobalID() const;

   inline RankType getLocalID() const; */

private:
  MPI_Comm gcomm_; // global communicator which includes manager

  MPI_Comm lcomm_;  // local communicator of pgroup

  RankType globalID_; // global rank

  RankType localID_; // rank inside ProcessGroup

  RankType managerID_; // rank of manager

  RankType localRootID_; // local root of pgroup

  TaskContainer tasks_; // task storage

  Task* currentTask_;

  StatusType status_;

  FullGrid<complex>* combinedFG_;

  DistributedSparseGridUniform<CombiDataType>* combinedUniDSG_;

  bool combinedFGexists_;

  CombiParameters combiParameters_;

  bool combiParametersSet_;

  void setCombinedSolutionUniform( Task* t );
};

inline Task* ProcessGroupWorker::getCurrentTask() {
  return currentTask_;
}

} /* namespace combigrid */

#endif /* PROCESSGROUPWORKER_HPP_ */

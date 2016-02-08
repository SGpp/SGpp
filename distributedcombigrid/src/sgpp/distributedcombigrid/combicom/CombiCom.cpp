/*
 * CombiCom.cpp
 *
 *  Created on: Oct 20, 2015
 *      Author: heenemo
 */
#include "sgpp/distributedcombigrid/combicom/CombiCom.hpp"

namespace combigrid {

static MPI_Comm distributedGlobalReduceCommunicator_ = MPI_COMM_NULL;

void CombiCom::registerDistributedGlobalReduceCommmunicator(int pgroupSize) {
  if (distributedGlobalReduceCommunicator_ != MPI_COMM_NULL)
    return;

  int globalID, globalSize;
  MPI_Comm_rank( MPI_COMM_WORLD, &globalID);
  MPI_Comm_size( MPI_COMM_WORLD, &globalSize);
  const int managerIDworld = globalSize - 1;

  // create communicator which only contains workers
  MPI_Comm workerComm;
  {
    int color = (globalID != managerIDworld) ? 1 : 0;
    int key = (globalID != managerIDworld) ? globalID : 0;
    MPI_Comm_split( MPI_COMM_WORLD, color, key, &workerComm);
  }

  if (globalID != managerIDworld) {
    int workerID;
    MPI_Comm_rank(workerComm, &workerID);

    MPI_Comm globalReduceComm;
    int color = workerID % pgroupSize;
    int key = workerID / pgroupSize;
    MPI_Comm_split(workerComm, color, key, &globalReduceComm);

    distributedGlobalReduceCommunicator_ = globalReduceComm;
  }
}

MPI_Comm CombiCom::getDistributedGlobalReduceCommunicator() {
  return distributedGlobalReduceCommunicator_;
}

} // namespace

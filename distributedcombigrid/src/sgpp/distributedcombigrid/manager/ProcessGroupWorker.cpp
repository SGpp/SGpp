/*
 * ProcessGroupWorker.cpp
 *
 *  Created on: Jun 24, 2014
 *      Author: heenemo
 */

#include "sgpp/distributedcombigrid/manager/ProcessGroupWorker.hpp"

#include "boost/lexical_cast.hpp"

#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupSignals.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGrid.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGridUniform.hpp"
#include "sgpp/distributedcombigrid/combicom/CombiCom.hpp"
#include "sgpp/distributedcombigrid/hierarchization/DistributedHierarchization.hpp"
#include "sgpp/distributedcombigrid/mpi/MPIUtils.hpp"

namespace combigrid {

ProcessGroupWorker::ProcessGroupWorker(MPI_Comm gcomm, MPI_Comm lcomm,
                                       RankType grank, RankType lrank, RankType manager, RankType lroot) :
  gcomm_(gcomm), lcomm_(lcomm), globalID_(grank), localID_(lrank), managerID_(
    manager), localRootID_(lroot), currentTask_( NULL), status_(
      PROCESS_GROUP_WAIT), combinedFG_( NULL), combinedFGexists_(false),
  combiParametersSet_(
    false), combinedUniDSG_( NULL) {
}

ProcessGroupWorker::~ProcessGroupWorker() {
  // TODO Auto-generated destructor stub
  delete combinedFG_;
}

SignalType ProcessGroupWorker::wait() {
  //std::cout << "rank " << globalID_ << " entering wait" << std::endl;

  if (status_ != PROCESS_GROUP_WAIT)
    return RUN_NEXT;

  SignalType signal = -1;

  if (localID_ == localRootID_) {
    // receive signal from manager
    MPI_Recv(&signal, 1, MPI_INT, managerID_, signalTag, gcomm_,
             MPI_STATUS_IGNORE);
  }

  // distribute signal to other processes of pgroup
  MPI_Bcast(&signal, 1, MPI_INT, localRootID_, lcomm_);

  std::cout << "received signal " << signal << std::endl;

  // process signal
  if (signal == RUN_FIRST) {
    Task* t;

    // local root receives task
    if (localID_ == localRootID_) {
      Task::receive(&t, managerID_, gcomm_);
    }

    // broadcast task to other process of pgroup
    Task::broadcast(&t, localRootID_, lcomm_);

    MPI_Barrier(lcomm_);

    // add task to task storage
    tasks_.push_back(t);

    status_ = PROCESS_GROUP_BUSY;

    // set currentTask
    currentTask_ = tasks_.back();

    // initalize task
    currentTask_->init(lcomm_);

    // execute task
    currentTask_->run(lcomm_);
  } else if (signal == RUN_NEXT) {
    // this should not happen
    assert(tasks_.size() > 0);

    // reset finished status of all tasks
    for (size_t i = 0; i < tasks_.size(); ++i)
      tasks_[i]->setFinished(false);

    status_ = PROCESS_GROUP_BUSY;

    // set currentTask
    currentTask_ = tasks_[0];

    // run first task
    currentTask_->run(lcomm_);
  } else if (signal == ADD_TASK) {
    std::cout << "adding a single task" << std::endl;

    Task* t;

    // local root receives task
    if (localID_ == localRootID_) {
      Task::receive(&t, managerID_, gcomm_);
    }

    // broadcast task to other process of pgroup
    Task::broadcast(&t, localRootID_, lcomm_);

    MPI_Barrier(lcomm_);

    /*
    if (localID_ == localRootID_)
      std::cout << "received new task with id = " << t->getID()
                << " and l = " << t->getLevelVector() << std::endl;
    */

    // check if task already exists on this group
    for ( auto tmp : tasks_ )
      assert( tmp->getID() != t->getID() );

    // initalize task
    t->init(lcomm_);

    t->setFinished( true );

    // add task to task storage
    tasks_.push_back(t);

    status_ = PROCESS_GROUP_BUSY;
  } else if (signal == EVAL) {
    // receive x

    // loop over all tasks
    // t.eval(x)
  } else if (signal == EXIT) {
    std::cout << "received exit signal" << std::endl;
  } else if (signal == SYNC_TASKS) {
    if (localID_ == localRootID_) {
      for (size_t i = 0; i < tasks_.size(); ++i) {
        Task::send(&tasks_[i], managerID_, gcomm_);
      }
    }
  } else if (signal == COMBINE) {
    std::cout << "starting distributed combination" << std::endl;
    combineUniform();
  } else if (signal == GRID_EVAL) {
    std::cout << "starting grid eval" << std::endl;
    gridEval();
  } else if (signal == COMBINE_FG) {
    std::cout << "starting fg combination" << std::endl;
    combineFG();
  } else if (signal == UPDATE_COMBI_PARAMETERS) {
    std::cout << "updating combi parameters " << std::endl;
    updateCombiParameters();
  } else if (signal == RECOMPUTE) {
    Task* t;

    // local root receives task
    if (localID_ == localRootID_) {
      Task::receive(&t, managerID_, gcomm_);
    }

    // broadcast task to other process of pgroup
    Task::broadcast(&t, localRootID_, lcomm_);

    MPI_Barrier(lcomm_);

    // add task to task storage
    tasks_.push_back(t);

    status_ = PROCESS_GROUP_BUSY;

    // set currentTask
    currentTask_ = tasks_.back();

    // initalize task
    currentTask_->init(lcomm_);

    // fill task with combisolution
    setCombinedSolutionUniform( currentTask_ );

    // execute task
    currentTask_->run(lcomm_);
  }


  // in the general case: send ready signal.
  ready();

  //std::cout << "leaving wait" << std::endl;

  return signal;
}

void ProcessGroupWorker::ready() {
  // check if there are unfinished tasks
  for (size_t i = 0; i < tasks_.size(); ++i) {
    if (!tasks_[i]->isFinished()) {
      status_ = PROCESS_GROUP_BUSY;

      // set currentTask
      currentTask_ = tasks_[i];
      currentTask_->run(lcomm_);
    }
  }

  // all tasks finished -> group waiting
  status_ = PROCESS_GROUP_WAIT;

  // send ready status to manager
  if (localID_ == localRootID_) {
    StatusType status = status_;
    MPI_Send(&status, 1, MPI_INT, managerID_, statusTag, gcomm_);
  }

  // reset current task
  currentTask_ = NULL;
}

void ProcessGroupWorker::combine() {
  // early exit if no tasks available
  // todo: doesnt work, each pgrouproot must call reduce function
  assert(tasks_.size() > 0);

  assert( combiParametersSet_ );
  DimType dim = combiParameters_.getDim();
  const LevelVector& lmin = combiParameters_.getLMin();
  const LevelVector& lmax = combiParameters_.getLMax();
  const std::vector<bool>& boundary = combiParameters_.getBoundary();

  // erzeug dsg
  DistributedSparseGrid<CombiDataType> dsg(dim, lmax, lmin, boundary, lcomm_);

  for (Task* t : tasks_) {
    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();

    // hierarchize dfg
    DistributedHierarchization::hierarchize<CombiDataType>(dfg);

    // lokales reduce auf sg ->
    //CombiCom::distributedLocalReduce<CombiDataType>( dfg, dsg, combiParameters_.getCoeff( t->getID() ) );
  }

  // globales reduce
  CombiCom::distributedGlobalReduce(dsg);

  for (Task* t : tasks_) {
    // get handle to dfg
    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();

    // lokales scatter von dsg auf dfg
    //CombiCom::distributedLocalScatter<CombiDataType>( dfg, dsg );

    // dehierarchize dfg
    DistributedHierarchization::dehierarchize<CombiDataType>(dfg);
  }
}

void ProcessGroupWorker::combineUniform() {
  // early exit if no tasks available
  // todo: doesnt work, each pgrouproot must call reduce function
  assert(tasks_.size() > 0);

  assert( combiParametersSet_ );
  DimType dim = combiParameters_.getDim();
  const LevelVector& lmin = combiParameters_.getLMin();
  LevelVector lmax = combiParameters_.getLMax();
  const std::vector<bool>& boundary = combiParameters_.getBoundary();;

  for (size_t i = 0; i < lmax.size(); ++i)
    if (lmax[i] > 1)
      lmax[i] -= 1;

  // todo: delete old dsg
  if (combinedUniDSG_ != NULL)
    delete combinedUniDSG_;

  theStatsContainer()->setTimerStart("combine_create_dsg");
  // erzeug dsg
  combinedUniDSG_ = new DistributedSparseGridUniform<CombiDataType>(dim, lmax,
      lmin, boundary,
      lcomm_);
  theStatsContainer()->setTimerStop("combine_create_dsg");

  // todo: move to init function to avoid reregistering
  // register dsg in all dfgs
  for (Task* t : tasks_) {
    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();

    dfg.registerUniformSG(*combinedUniDSG_);
  }

  real thierarchization = 0.0;
  real taddtosg = 0.0;

  theStatsContainer()->setTimerStart("combine_local_reduce");

  for (Task* t : tasks_) {
    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();

    real thierarchization_start = MPI_Wtime();
    // hierarchize dfg
    DistributedHierarchization::hierarchize<CombiDataType>(dfg);
    thierarchization += MPI_Wtime() - thierarchization_start;

    real taddtosg_start = MPI_Wtime();
    // lokales reduce auf sg ->
    dfg.addToUniformSG( *combinedUniDSG_, combiParameters_.getCoeff( t->getID() ) );
    taddtosg += MPI_Wtime() - taddtosg_start;
  }

  theStatsContainer()->setTimerStop("combine_local_reduce");

  theStatsContainer()->setValue("combine_local_reduce_hierarchization",
                                thierarchization);
  theStatsContainer()->setValue("combine_local_reduce_addtosg", taddtosg);

  theStatsContainer()->setTimerStart("combine_global_reduce");
  CombiCom::distributedGlobalReduce( *combinedUniDSG_ );
  theStatsContainer()->setTimerStop("combine_global_reduce");

  theStatsContainer()->setTimerStart("combine_local_scatter");

  for (Task* t : tasks_) {
    // get handle to dfg
    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();

    // extract dfg vom dsg
    dfg.extractFromUniformSG( *combinedUniDSG_ );

    // dehierarchize dfg
    DistributedHierarchization::dehierarchize<CombiDataType>( dfg );
  }

  theStatsContainer()->setTimerStop("combine_local_scatter");
}

void ProcessGroupWorker::gridEval() {
  /* error if no tasks available
   * todo: however, this is not a real problem, we could can create an empty
   * grid an contribute to the reduce operation. at the moment even the dim
   * parameter is stored in the tasks, so if no task available we have no access
   * to this parameter.
   */
  assert(tasks_.size() > 0);

  assert(combiParametersSet_);
  const DimType dim = combiParameters_.getDim();

  LevelVector leval(dim);

  // receive leval
  if (localID_ == localRootID_) {
    // receive size of levelvector = dimensionality
    MPI_Status status;
    int bsize;
    MPI_Probe(managerID_, 0, gcomm_, &status);
    MPI_Get_count(&status, MPI_INT, &bsize);

    assert(bsize == static_cast<int>(dim));

    std::vector<int> tmp(dim);
    MPI_Recv(&tmp[0], bsize, MPI_INT, managerID_, 0, gcomm_, MPI_STATUS_IGNORE);
    leval = LevelVector(tmp.begin(), tmp.end());
  }

  assert( combiParametersSet_ );
  const std::vector<bool>& boundary = combiParameters_.getBoundary();
  FullGrid<CombiDataType> fg_red(dim, leval, boundary);

  // create the empty grid on only on localroot
  if (localID_ == localRootID_) {
    fg_red.createFullGrid();
  }

  // collect fg on pgrouproot and reduce
  for (size_t i = 0; i < tasks_.size(); ++i) {
    Task* t = tasks_[i];

    FullGrid<CombiDataType> fg(t->getDim(), t->getLevelVector(), boundary );

    if (localID_ == localRootID_) {
      fg.createFullGrid();
    }

    t->getFullGrid(fg, localRootID_, lcomm_);

    // local root section
    if (localID_ == localRootID_) {
      fg_red.add(fg, combiParameters_.getCoeff( t->getID() ) );
    }
  }

  // global reduce of f_red
  if (localID_ == localRootID_)
    CombiCom::FGReduce(fg_red, managerID_, gcomm_);
}

//todo: this is just a temporary function which will drop out some day
// also this function requires a modified fgreduce method which uses allreduce
// instead reduce in manger
void ProcessGroupWorker::combineFG() {
  //gridEval();

  // TODO: Sync back to fullgrids
}

void ProcessGroupWorker::updateCombiParameters() {
  CombiParameters tmp;

  // local root receives task
  if (localID_ == localRootID_) {
    MPIUtils::receiveClass(&tmp, managerID_, gcomm_);
  }

  // broadcast task to other process of pgroup
  MPIUtils::broadcastClass(&tmp, localRootID_, lcomm_);

  combiParameters_ = tmp;

  combiParametersSet_ = true;
}


void ProcessGroupWorker::setCombinedSolutionUniform( Task* t ) {
  assert( combinedUniDSG_ != NULL );

  // get handle to dfg
  DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();

  // extract dfg vom dsg
  dfg.extractFromUniformSG( *combinedUniDSG_ );

  // dehierarchize dfg
  DistributedHierarchization::dehierarchize<CombiDataType>( dfg );
}

} /* namespace combigrid */

/*
 * combi_example.cpp
 *
 *  Created on: Sep 23, 2015
 *      Author: heenemo
 */
#include <mpi.h>
#include <vector>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/serialization/export.hpp>

// compulsory includes for basic functionality
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupWorker.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessManager.hpp"

// include user specific task. this is the interface to your application
#include "TaskExample.hpp"

using namespace combigrid;

// this is necessary for correct function of task serialization
BOOST_CLASS_EXPORT(TaskExample)

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  // read in parameter file
  boost::property_tree::ptree cfg;
  boost::property_tree::ini_parser::read_ini("ctparam", cfg);

  // number of process groups and number of processes per group
  size_t ngroup = cfg.get<size_t>("manager.ngroup");
  size_t nprocs = cfg.get<size_t>("manager.nprocs");

  //TODO: put creation of communicators in distinct functions

  /* determine global rank of each process
   * the manager process always has the highest rank
   * all other processes are worker processes */
  int globalID, globalSize;
  MPI_Comm_rank(MPI_COMM_WORLD, &globalID);
  MPI_Comm_size(MPI_COMM_WORLD, &globalSize);
  const int managerIDworld = globalSize - 1;

  assert(globalSize == int(ngroup * nprocs + 1));

  /* create a local communicator for each process group
   * lcomm is the local communicator of its own process group for each worker process
   * for manager, lcomm is a group which contains only manager process and can be ignored
   */
  MPI_Comm lcomm;
  int color = globalID / nprocs;
  int key = globalID - color * nprocs;
  int lrank;
  MPI_Comm_split(MPI_COMM_WORLD, color, key, &lcomm);
  MPI_Comm_rank(lcomm, &lrank);

  /* create global communicator which contains only the manager and the master
   * process of each process group
   * the master processes of the process groups are the processes which have
   * rank 0 in lcomm
   * this communicator is used for communication between master processes of the
   * process groups and the manager and the master processes to each other
   */
  MPI_Group worldGroup;
  MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);

  std::vector<int> ranks(ngroup + 1);

  for (size_t i = 0; i < ngroup; i++) {
    ranks[i] = i * nprocs;
  }

  ranks.back() = managerIDworld;

  MPI_Group rootGroup;
  MPI_Group_incl(worldGroup, (int) ranks.size(), &ranks[0], &rootGroup);

  MPI_Comm gcomm;
  MPI_Comm_create(MPI_COMM_WORLD, rootGroup, &gcomm);

  int managerIDgcomm = MPI_PROC_NULL;
  int grank;

  if (gcomm != MPI_COMM_NULL) {
    int gcommSize;
    MPI_Comm_size(gcomm, &gcommSize);
    managerIDgcomm = gcommSize - 1;
    MPI_Comm_rank(gcomm, &grank);
  }

  CombiCom::registerDistributedGlobalReduceCommmunicator(nprocs);

  // manager code
  if (globalID == managerIDworld) {
    /* create an abstraction of the process groups for the manager's view
     * a pgroup is identified by the ID in gcomm
     */
    ProcessGroupManagerContainer pgroups;

    for (size_t i = 0; i < ngroup; ++i) {
      int pgroupRootID(i);
      ProcessGroupManager grp(managerIDgcomm, pgroupRootID, gcomm);
      pgroups.push_back(grp);
    }

    // create load model
    LoadModel* loadmodel = new LinearLoadModel();

    /* generate a list of levelvectors and coefficients
     * CombiTS_CT will generate a valid combination. however, you could
     * also read in a list of levelvectors and coefficients from a file */
    DimType dim = cfg.get<DimType>("ct.dim");
    LevelVector lmin(dim), lmax(dim), leval(dim);
    IndexVector p(dim);
    real dt;
    size_t nsteps, ncombi;
    cfg.get<std::string>("ct.lmin") >> lmin;
    cfg.get<std::string>("ct.lmax") >> lmax;
    cfg.get<std::string>("ct.leval") >> leval;
    cfg.get<std::string>("ct.p") >> p;
    ncombi = cfg.get<size_t>("ct.ncombi");
    dt = cfg.get<real>("application.dt");
    nsteps = cfg.get<size_t>("application.nsteps");

    // todo: read in boundary vector from ctparam
    std::vector<bool> boundary(dim, true);

    // check parallelization vector p agrees with nprocs
    IndexType checkProcs = 1;

    for (auto k : p)
      checkProcs *= k;

    assert(checkProcs == IndexType(nprocs));

    CombiMinMaxScheme combischeme(dim, lmin, lmax);
    std::vector<LevelVector> levels = combischeme.getCombiSpaces();
    std::vector<real> coeffs = combischeme.getCoeffs();

    // output of combination setup
    std::cout << "lmin = " << lmin << std::endl;
    std::cout << "lmax = " << lmax << std::endl;
    std::cout << "CombiScheme: " << std::endl;
    std::cout << combischeme << std::endl;



    // create Tasks
    TaskContainer tasks;
    std::vector<int> taskIDs;

    for (size_t i = 0; i < levels.size(); i++) {
      Task* t = new TaskExample(dim, levels[i], lmax, lmin, boundary, coeffs[i],
                                loadmodel, dt, nsteps, p);
      tasks.push_back(t);
      taskIDs.push_back( t->getID() );
    }

    // create combiparamters
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs );

    // create Manager with process groups
    ProcessManager manager(pgroups, tasks, params, managerIDgcomm, gcomm);

    // combiparameters need to be set before starting the computation
    manager.updateCombiParameters();

    /* distribute task according to load model and start computation for
     * the first time */
    manager.runfirst();

    std::ofstream myfile;
    myfile.open("out/solution.dat");

    for (size_t i = 0; i < ncombi; ++i) {
      manager.combine();

      // evaluate solution
      FullGrid<double> fg_eval(dim, leval, boundary);
      manager.gridEval(fg_eval);

      // write solution to file
      std::vector<double> coords(dim, 0.0);

      for (int i = 0; i < fg_eval.getNrElements(); i++) {
        if (i % fg_eval.length(0) == 0 && i > 0) {
          myfile << std::endl;
        }

        fg_eval.getCoords(i, coords);
        myfile << coords[0] << "\t" << coords[1] << "\t"
               << fg_eval.getElementVector()[i] << std::endl;
      }

      myfile << std::endl << std::endl;

      // run tasks for next time interval
      manager.runnext();
    }

    myfile.close();

    // send exit signal to workers in order to enable a clean program termination
    manager.exit();
  }

  // worker code
  else {
    // create abstraction of the process group from the worker's view
    ProcessGroupWorker pgroup(gcomm, lcomm, grank, lrank, managerIDgcomm);

    // wait for instructions from manager
    SignalType signal = -1;

    while (signal != EXIT)
      signal = pgroup.wait();
  }

  MPI_Finalize();

  return 0;
}

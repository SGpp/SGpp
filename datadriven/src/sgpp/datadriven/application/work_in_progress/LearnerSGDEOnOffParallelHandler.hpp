// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_LEARNERSGDEONOFFPARALLELHANDLER_H
#define SGPP_LEARNERSGDEONOFFPARALLELHANDLER_H

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/algorithm/ConvergenceMonitor.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/application/LearnerSGDEOnOff.hpp>
#include <sgpp/datadriven/application/work_in_progress/MPITaskScheduler.hpp>
#include <sgpp/datadriven/application/work_in_progress/LearnerSGDEOnOffParallel.hpp>

#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <mpi.h>


namespace sgpp {
namespace datadriven {

class LearnerSGDEOnOffParallelHandler {
 protected:
  std::vector<RefinementResult> vectorRefinementResults;
  LearnerSGDEOnOffParallel *learnerInstance;

  unsigned long
  handleDataAndZeroBasedRefinement(bool preCompute,
                                   MultiGridRefinementFunctor *func,
                                   unsigned long idx,
                                   base::Grid &grid,
                                   base::GridGenerator &gridGen) const;

  unsigned long
  handleSurplusBasedRefinement(DBMatOnlineDE *densEst,
                               Grid &grid,
                               base::GridGenerator &gridGen) const;

 public:
  LearnerSGDEOnOffParallelHandler(LearnerSGDEOnOffParallel *learnerInstance,
                                  size_t numClasses);

  /**
      * After refinement completes locally or refinement results have been received over MPI, this method uses the results to adjust the grid and alpha vector.
      * If this is run on the master, the grid changes will be exported over MPI.
      * If the system matrix is refineable, a system matrix update will be assigned to workers over MPI.
      * @param classIndex The index of the class being updated
      * @param refinementResult The grid changes from the refinement cycle
      * @param densEst A pointer to the online object specfic to this class
      */
  void updateClassVariablesAfterRefinement(unsigned long classIndex,
                                           RefinementResult *refinementResult,
                                           DBMatOnlineDE *densEst);

  /**
       * Fetches the currently stored refinement results for a specific class
       * @param classIndex The class to search refinement results for
       * @return A reference to the stored refinement results
       */
  RefinementResult &getRefinementResult(unsigned long classIndex);

/**
    * Check whether all grids are consistent and the scheduler is currently allowing refinement.
    * @return Whether refinement is currently possible
    */
  bool checkReadyForRefinement() const;

  bool checkRefinementNecessary(const std::string &refMonitor, unsigned long refPeriod,
                                unsigned long totalInstances,
                                double currentValidError, double currentTrainError,
                                unsigned long numberOfCompletedRefinements,
                                ConvergenceMonitor &monitor);

  /**
           * Handles refinement for a specific class.
           * @param refType
           * @param refinementResult
           * @param onlineObjects
           * @param preCompute
           * @param refinementFunctor
           * @param classIndex
           */
  void doRefinementForClass(const std::string &refType,
                            RefinementResult *refinementResult,
                            const ClassDensityConntainer &onlineObjects,
                            bool preCompute,
                            MultiGridRefinementFunctor *refinementFunctor,
                            unsigned long classIndex);
};
}  // namespace datadriven
}  // namespace sgpp

#endif  // SGPP_LEARNERSGDEONOFFPARALLELHANDLER_H

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/AuxiliaryStructures.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/MPITaskScheduler.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitor.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorPeriodic.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorConvergence.hpp>

#include <mpi.h>

#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace sgpp {
namespace datadriven {

// Forward declare learner, as we use only pointer
class LearnerSGDEOnOffParallel;

class RefinementHandler {
 protected:
  std::vector<RefinementResult> vectorRefinementResults;
  LearnerSGDEOnOffParallel *learnerInstance;

  /**
   * Logic that handles data-based and zero-crossing refinement functors
   * @param preCompute Whether to precompute evaluations in the functor
   * @param func Pointer to the refinement functor itself
   * @param idx Class index
   * @param grid The grid for the current class
   * @param gridGen The grid's generator for the current grid
   * @return The number of added grid points
   */
  size_t handleDataAndZeroBasedRefinement(bool preCompute, MultiGridRefinementFunctor *func,
                                          size_t idx, base::Grid &grid,
                                          base::GridGenerator &gridGen) const;

  /**
   * Logic that handles surplus based refinement functors
   * @param densEst Online objects for use in density estimation for the current class
   * @param grid The current classes grid
   * @param alpha The current surpluss vector
   * @param gridGen The current grid's grid generator
   * @param adaptivityConfig the configuration for the adaptivity
   * @return The number of added grid points
   */
  size_t handleSurplusBasedRefinement(DBMatOnlineDE *densEst, Grid &grid, DataVector &alpha,
                                      base::GridGenerator &gridGen,
                                      sgpp::base::AdaptivityConfiguration adaptivityConfig) const;

 public:
  /**
   * Creates the refinement handler for a specific learner instance
   * @param learnerInstance The instance of the learner to handle refinement for
   * @param numClasses The number of classes for the current problem
   */
  RefinementHandler(LearnerSGDEOnOffParallel *learnerInstance, size_t numClasses);

  /**
    * After refinement completes locally or refinement results have been received over MPI, this
   * method uses the results to adjust the grid and alpha vector.
    * If this is run on the master, the grid changes will be exported over MPI.
    * If the system matrix is refineable, a system matrix update will be assigned to workers over
   * MPI.
    * @param classIndex The index of the class being updated
    * @param refinementResult The grid changes from the refinement cycle
    * @param densEst A pointer to the online object specfic to this class
    * @param grid the grid of the current density function
    */
  void updateClassVariablesAfterRefinement(size_t classIndex, RefinementResult *refinementResult,
                                           DBMatOnlineDE *densEst, Grid &grid);

  /**
   * Fetches the currently stored refinement results for a specific class
   * @param classIndex The class to search refinement results for
   * @return A reference to the stored refinement results
   */
  RefinementResult &getRefinementResult(size_t classIndex);

  /**
    * Check whether all grids are consistent and the scheduler is currently allowing refinement.
    * @return Whether refinement is currently possible
    */
  bool checkReadyForRefinement() const;

  /**
   * Check whether refinement is currently necessary according to the guidelines set by the user
   * @param refMonitor String constant specifying the monitor to use for refinement
   * @param refPeriod The minimum period in which refinement cycles are allowed
   * @param batchSize The number of instances that were added by the current batch
   * @param currentValidError The current validation error
   * @param currentTrainError The current training error
   * @param numberOfCompletedRefinements The number of refinement cycles already completed
   * @param monitor The convergence monitor, if any
   * @param adaptivityConfig the configuration for the adaptivity of the grids
   * @return How many refinement cycles should be started
   */
  size_t checkRefinementNecessary(const std::string &refMonitor, size_t refPeriod, size_t batchSize,
                                  double currentValidError, double currentTrainError,
                                  size_t numberOfCompletedRefinements, RefinementMonitor &monitor,
                                  sgpp::base::AdaptivityConfiguration adaptivityConfig);

  /**
   * Handles refinement for a specific class.
   * @param refType String constant specifying the type of refinement functor
   * @param refinementResult The RefinementResult used to store changes for the grid
   * @param onlineObjects The density estimation online objects
   * @param grid the grid of the online object of the current class
   * @param alpha the surplusses of the current class
   * @param preCompute Whether to precompute the functor's evaluation step
   * @param refinementFunctor The refinement functor to use
   * @param classIndex The index of the current class for which refinement is taking place
   * @param adaptivityConfig the configuration for the adaptivity of the grids
   */
  void doRefinementForClass(
      const std::string &refType, RefinementResult *refinementResult,
      const std::vector<std::pair<std::unique_ptr<DBMatOnlineDE>, size_t>> &onlineObjects,
      Grid &grid, DataVector &alpha, bool preCompute, MultiGridRefinementFunctor *refinementFunctor,
      size_t classIndex, sgpp::base::AdaptivityConfiguration &adaptivityConfig);
};
}  // namespace datadriven
}  // namespace sgpp

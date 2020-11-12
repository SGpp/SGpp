// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp_base.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/LearnerSGDEOnOffParallel.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/MPIMethods.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/RefinementHandler.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitor.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorPeriodic.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorConvergence.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <vector>
#include <string>

namespace sgpp {
namespace datadriven {
bool RefinementHandler::checkReadyForRefinement() const {
  // All local grids in a consistent state and have OK from scheduler
  return learnerInstance->getScheduler().isReadyForRefinement() &&
         learnerInstance->checkAllGridsConsistent();
}

void RefinementHandler::doRefinementForClass(
    const std::string &refType, RefinementResult *refinementResult,
    const std::vector<std::pair<std::unique_ptr<DBMatOnlineDE>, size_t>> &onlineObjects, Grid &grid,
    DataVector &alpha, bool preCompute, MultiGridRefinementFunctor *refinementFunctor,
    size_t classIndex, sgpp::base::AdaptivityConfiguration &adaptivityConfig) {
  // perform refinement/coarsening for grid which corresponds to current
  // index
  std::cout << "Refinement and coarsening for class: " << classIndex << std::endl;
  auto densEst = onlineObjects[classIndex].first.get();

  size_t oldGridSize = grid.getSize();
  D(std::cout << "Size before adaptivity: " << oldGridSize << std::endl;)

  base::GridGenerator &gridGen = grid.getGenerator();

  size_t numberOfNewPoints = 0;

  if (refType == "surplus") {
    numberOfNewPoints =
        handleSurplusBasedRefinement(densEst, grid, alpha, gridGen, adaptivityConfig);
  } else if ((refType == "data") || (refType == "zero")) {
    numberOfNewPoints =
        handleDataAndZeroBasedRefinement(preCompute, refinementFunctor, classIndex, grid, gridGen);
  }

  size_t newGridSize = grid.getSize();
  std::cout << "grid size after adaptivity: " << newGridSize << " (previously " << oldGridSize
            << "), " << numberOfNewPoints << " new points on grid" << std::endl;

  if (numberOfNewPoints != newGridSize - oldGridSize) {
    std::cout << "Reported grid sizes do not match up (refined " << numberOfNewPoints << ", old "
              << oldGridSize << ", new " << newGridSize << ")" << std::endl;
    throw sgpp::base::algorithm_exception("Reported grid delta not valid.");
  }

  D(std::cout << "Preparing refinement result update" << std::endl;)
  D(std::cout << "Clearing old added grid points" << std::endl;)
  refinementResult->addedGridPoints.clear();

  D(std::cout << "Clearing old deleted grid points" << std::endl;)
  refinementResult->deletedGridPointsIndices.clear();

  size_t numDimensions = learnerInstance->getDimensionality();
  // Collect new grid points into the refinement result for shipping
  for (size_t i = 0; i < numberOfNewPoints; i++) {
    std::vector<LevelIndexPair> levelIndexVector(numDimensions);
    base::HashGridPoint &gridPoint = grid.getStorage()[oldGridSize + i];
    for (size_t currentDimension = 0; currentDimension < numDimensions; currentDimension += 1) {
      uint32_t pointLevel;
      uint32_t pointIndex;
      gridPoint.get(currentDimension, pointLevel, pointIndex);

      levelIndexVector[currentDimension].level = pointLevel;
      levelIndexVector[currentDimension].index = pointIndex;

      //                    std::cout << "Fetched level index vector: dimension " <<
      //                    currentDimension
      //                              << ", level " << pointLevel <<
      //                              ", index " << pointIndex << std::endl;
    }

    refinementResult->addedGridPoints.push_back(levelIndexVector);
  }

  updateClassVariablesAfterRefinement(classIndex, refinementResult, densEst, grid);
}

void RefinementHandler::updateClassVariablesAfterRefinement(size_t classIndex,
                                                            RefinementResult *refinementResult,
                                                            DBMatOnlineDE *densEst, Grid &grid) {
  if (!MPIMethods::isMaster()) {
    std::cout << "Applying refinement result class " << classIndex << " from master" << std::endl;
    D(std::cout << "Old grid " << classIndex << " size is " << grid.getSize() << std::endl;)

    size_t numDimensions = learnerInstance->getDimensionality();

    // Delete the grid points removed on master thread
    grid.getStorage().deletePoints(refinementResult->deletedGridPointsIndices);

    size_t sizeBeforeAdditions = grid.getSize();
    D(std::cout << "Grid size after deleting is " << sizeBeforeAdditions << std::endl;)

    // Add grid points added on master thread
    for (std::vector<LevelIndexPair> &levelIndexVector : refinementResult->addedGridPoints) {
      D(size_t sizeBeforePoint = grid.getSize();)
      auto *gridPoint = new sgpp::base::HashGridPoint(numDimensions);

      for (size_t currentDimension = 0; currentDimension < numDimensions; currentDimension++) {
        gridPoint->set(currentDimension, levelIndexVector[currentDimension].level,
                       levelIndexVector[currentDimension].index);
      }
      grid.getStorage().insert(*gridPoint);
      D(size_t sizeAfterPoint = grid.getSize(); if (sizeAfterPoint - sizeBeforePoint != 1) {
        std::cout << "Inserted grid point but size change incorrect (old " << sizeBeforePoint
                  << ", new " << sizeAfterPoint << "), point " << gridPoint->getHash() << std::endl;
        throw sgpp::base::algorithm_exception("Inserted grid point did not cause grid delta.");
      })
    }

    grid.getStorage().recalcLeafProperty();

    size_t sizeAfterAdditions = grid.getSize();
    D(std::cout << "New grid " << classIndex << " size is " << sizeAfterAdditions << std::endl;)
    if (sizeAfterAdditions - sizeBeforeAdditions != refinementResult->addedGridPoints.size()) {
      std::cout << "Grid growth not correlated to refinement results (grid delta "
                << sizeAfterAdditions - sizeBeforeAdditions << ", additions "
                << refinementResult->addedGridPoints.size() << ")" << std::endl;
      throw sgpp::base::algorithm_exception("Grid growth from refinement not correct.");
    }
  } else {
    // apply grid changes to the system matrix factorization
    size_t currentGridVersion = learnerInstance->getLocalGridVersion(classIndex);
    if (!learnerInstance->checkGridStateConsistent(classIndex)) {
      std::cout << "Attempting to increment grid version on non consistent grid " << classIndex
                << " version " << currentGridVersion << std::endl;
      throw sgpp::base::algorithm_exception("Setting of grid version failed due to inconsistent.");
    }

    learnerInstance->setLocalGridVersion(classIndex, currentGridVersion + 1);

    // Send class update in preparation for system matrix decomposition update
    MPIMethods::sendRefinementUpdates(classIndex, refinementResult->deletedGridPointsIndices,
                                      refinementResult->addedGridPoints);

    if (learnerInstance->getOffline()->isRefineable()) {
      AssignTaskResult assignTaskResult{};
      learnerInstance->getScheduler().assignTaskStaticTaskSize(
          RECOMPUTE_SYSTEM_MATRIX_DECOMPOSITION, assignTaskResult);

      std::cout << "Assigning update of system matrix decomposition " << classIndex << " to worker "
                << assignTaskResult.workerID << std::endl;
      MPIMethods::assignSystemMatrixUpdate(assignTaskResult.workerID, classIndex);
    }
  }

  // update alpha vector
  D(size_t oldSize = densEst->getAlpha().size();)
  learnerInstance->updateAlpha(classIndex, &(refinementResult->deletedGridPointsIndices),
                               refinementResult->addedGridPoints.size());
  D(std::cout << "Updated alpha vector " << classIndex << " (old size " << oldSize << ", new size "
              << densEst->getAlpha().size() << ")" << std::endl;)
}

size_t RefinementHandler::checkRefinementNecessary(
    const std::string &refMonitor, size_t refPeriod, size_t batchSize, double currentValidError,
    double currentTrainError, size_t numberOfCompletedRefinements, RefinementMonitor &monitor,
    sgpp::base::AdaptivityConfiguration adaptivityConfig) {
  auto &offline = learnerInstance->getOffline();
  // access DBMatOnlineDE-objects of all classes in order
  // to apply adaptivity to the specific sparse grids later on

  // check if and how many refinements should be performed
  size_t refinementsNecessary = 0;
  if (offline->isRefineable() && numberOfCompletedRefinements < adaptivityConfig.numRefinements_) {
    currentValidError = learnerInstance->getError(*learnerInstance->getValidationData());
    currentTrainError =
        learnerInstance->getError(learnerInstance->getTrainData());  // if train dataset is large
    // use a subset for error
    monitor.pushToBuffer(batchSize, currentValidError, currentTrainError);
    refinementsNecessary = monitor.refinementsNecessary();
  }
  return refinementsNecessary;
}

size_t RefinementHandler::handleDataAndZeroBasedRefinement(bool preCompute,
                                                           MultiGridRefinementFunctor *func,
                                                           size_t idx, base::Grid &grid,
                                                           base::GridGenerator &gridGen) const {
  if (preCompute) {
    // precompute the evals (needs to be done once per step, before
    // any refinement is done
    func->preComputeEvaluations();
  }
  func->setGridIndex(idx);
  // perform refinement (zero-crossings-based / data-based)
  size_t gridSizeBeforeRefine = grid.getSize();
  gridGen.refine(*func);
  size_t gridSizeAfterRefine = grid.getSize();
  return gridSizeAfterRefine - gridSizeBeforeRefine;
}

size_t RefinementHandler::handleSurplusBasedRefinement(
    DBMatOnlineDE *densEst, base::Grid &grid, DataVector &alpha, base::GridGenerator &gridGen,
    sgpp::base::AdaptivityConfiguration adaptivityConfig) const {
  DataVector *alphaWork;  // required for surplus refinement
  // auxiliary variables
  DataVector p(learnerInstance->getTrainData().getDimension());

  std::unique_ptr<sgpp::base::OperationEval> opEval(op_factory::createOperationEval(grid));
  sgpp::base::HashGridStorage &gridStorage = grid.getStorage();
  alphaWork = &alpha;
  DataVector alphaWeight(alphaWork->getSize());
  // determine surpluses
  for (size_t k = 0; k < gridStorage.getSize(); k++) {
    // sets values of p to the coordinates of the given GridPoint gp
    gridStorage.getPoint(k).getStandardCoordinates(p);
    // multiply k-th alpha with the evaluated function at grind-point
    // k
    alphaWeight[k] = alphaWork->get(k) * opEval->eval(*alphaWork, p);
  }

  // Perform Coarsening (surplus based)
  /*if (coarseCnt < maxCoarseNum) {
    HashCoarsening coarse_;
    // std::cout << "\n" << "Start coarsening\n";

    // Coarsening based on surpluses
    SurplusCoarseningFunctor scf(
      alphaWeight, coarseNumPoints, coarseThreshold);

    // std::cout << "Size before coarsening:" << grid->getSize() <<
  "\n";
    // int old_size = grid->getSize();
    coarse_.free_coarsen_NFirstOnly(
      grid->getStorage(), scf, alphaWeight, grid->getSize());

    std::cout << "Size after coarsening:" << grid->getSize() <<
  "\n\n";
    // int new_size = grid->getSize();

    deletedGridPoints.clear();
    deletedGridPoints = coarse_.getDeletedPoints();

    (*refineCoarse)[idx].first = deletedGridPoints;

    coarseCnt++;
  }*/

  // perform refinement (surplus based)
  size_t sizeBeforeRefine = grid.getSize();
  // simple refinement based on surpluses
  sgpp::base::SurplusRefinementFunctor srf(alphaWeight, adaptivityConfig.numRefinementPoints_);
  gridGen.refine(srf);
  size_t sizeAfterRefine = grid.getSize();
  return sizeAfterRefine - sizeBeforeRefine;
}

RefinementResult &RefinementHandler::getRefinementResult(size_t classIndex) {
  return vectorRefinementResults[classIndex];
}

RefinementHandler::RefinementHandler(LearnerSGDEOnOffParallel *learnerInstance, size_t numClasses) {
  this->learnerInstance = learnerInstance;
  RefinementResult initResult{};
  vectorRefinementResults.insert(vectorRefinementResults.begin(), numClasses, initResult);
}
}  // namespace datadriven
}  // namespace sgpp

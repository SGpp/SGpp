//
// Created by Vincent_Bode on 03.09.2017.
//

#include <sgpp_base.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/algorithm/ConvergenceMonitor.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/application/work_in_progress/MPIMethods.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>


namespace sgpp {
namespace datadriven {
bool LearnerSGDEOnOffParallelHandler::checkReadyForRefinement() const {
    // All local grids in a consistent state and have OK from scheduler
    return learnerInstance->getScheduler().isReadyForRefinement()
           && learnerInstance->checkAllGridsConsistent();
}

void LearnerSGDEOnOffParallelHandler::doRefinementForClass(const std::string &refType,
                                                           RefinementResult *refinementResult,
                                                           const ClassDensityConntainer &onlineObjects,
                                                           bool preCompute,
                                                           MultiGridRefinementFunctor *refinementFunctor,
                                                           unsigned long classIndex) {
    // perform refinement/coarsening for grid which corresponds to current
    // index
    std::cout << "Refinement and coarsening for class: " << classIndex
              << std::endl;
    auto densEst = onlineObjects[classIndex].first.get();
    Grid &grid = densEst->getOfflineObject().getGrid();

    unsigned long oldGridSize = grid.getSize();
    D(std::cout << "Size before adaptivity: " << oldGridSize
                << std::endl;)

    base::GridGenerator &gridGen = grid.getGenerator();

    unsigned long numberOfNewPoints = 0;

    if (refType == "surplus") {
        numberOfNewPoints = handleSurplusBasedRefinement(densEst, grid, gridGen);

    } else if ((refType == "data") || (refType == "zero")) {
        numberOfNewPoints = handleDataAndZeroBasedRefinement(preCompute, refinementFunctor,
                                                             classIndex, grid,
                                                             gridGen);
    }

    unsigned long newGridSize = grid.getSize();
    std::cout << "grid size after adaptivity: " << newGridSize
              << " (previously " << oldGridSize
              << "), " << numberOfNewPoints << " new points on grid"
              << std::endl;

    if (numberOfNewPoints != newGridSize - oldGridSize) {
        std::cout << "Reported grid sizes do not match up (refined " << numberOfNewPoints
                  << ", old "
                  << oldGridSize << ", new " << newGridSize << ")" << std::endl;
        exit(-1);
    }

    D(std::cout << "Preparing refinement result update" << std::endl;)
    D(std::cout << "Clearing old added grid points" << std::endl;)
    refinementResult->addedGridPoints.clear();

    D(std::cout << "Clearing old deleted grid points" << std::endl;)
    refinementResult->deletedGridPointsIndexes.clear();


    unsigned long numDimensions = learnerInstance->getDimensionality();
    //Collect new grid points into the refinement result for shipping
    for (unsigned long i = 0; i < numberOfNewPoints; i++) {
        std::vector<LevelIndexPair> levelIndexVector(numDimensions);
        base::HashGridPoint &gridPoint = grid.getStorage()[oldGridSize + i];
        for (unsigned long currentDimension = 0;
             currentDimension < numDimensions; currentDimension += 1) {
            uint32_t pointLevel;
            uint32_t pointIndex;
            gridPoint.get(currentDimension, pointLevel, pointIndex);

            levelIndexVector[currentDimension].level = pointLevel;
            levelIndexVector[currentDimension].index = pointIndex;

//                    std::cout << "Fetched level index vector: dimension " << currentDimension
//                              << ", level " << pointLevel <<
//                              ", index " << pointIndex << std::endl;
        }

        refinementResult->addedGridPoints.push_back(levelIndexVector);

        D(printPoint(gridPoint))
    }

    updateClassVariablesAfterRefinement(classIndex, refinementResult, densEst);

}

void LearnerSGDEOnOffParallelHandler::updateClassVariablesAfterRefinement(unsigned long classIndex,
                                                                          RefinementResult *refinementResult,
                                                                          DBMatOnlineDE *densEst) {

    base::Grid &grid = densEst->getOfflineObject().getGrid();

    if (!MPIMethods::isMaster()) {
        std::cout << "Applying refinement result class " << classIndex << " from master"
                  << std::endl;
        D(std::cout << "Old grid " << classIndex << " size is " << grid.getSize() << std::endl;)


        unsigned long numDimensions = learnerInstance->getDimensionality();

        //Delete the grid points removed on master thread
        grid.getStorage().deletePoints(refinementResult->deletedGridPointsIndexes);

        unsigned long sizeBeforeAdditions = grid.getSize();
        D(std::cout << "Grid size after deleting is " << sizeBeforeAdditions << std::endl;)

        //Add grid points added on master thread
        for (std::vector<LevelIndexPair> &levelIndexVector : refinementResult->addedGridPoints) {
            D(unsigned long sizeBeforePoint = grid.getSize();)
            auto *gridPoint = new sgpp::base::HashGridPoint(numDimensions);
            //TODO: What happens when other points are changed (ie Leaf boolean etc)

            for (unsigned long currentDimension = 0;
                 currentDimension < numDimensions; currentDimension++) {
                gridPoint->set(currentDimension,
                               levelIndexVector[currentDimension].level,
                               levelIndexVector[currentDimension].index);
            }
            grid.getStorage().insert(*gridPoint);
            D(
                    unsigned long
                            sizeAfterPoint = grid.getSize();
                    if (sizeAfterPoint - sizeBeforePoint != 1) {
                        std::cout << "Inserted grid point but size change incorrect (old "
                                  << sizeBeforePoint
                                  << ", new " << sizeAfterPoint << "), point "
                                  << gridPoint->getHash()
                                  << std::endl;
                        printPoint(gridPoint);
                        exit(-1);
                    }
            )
        }

        //TODO: This might be unnecessary
        grid.getStorage().recalcLeafProperty();

        unsigned long sizeAfterAdditions = grid.getSize();
        D(std::cout << "New grid " << classIndex << " size is " << sizeAfterAdditions << std::endl;)
        if (sizeAfterAdditions - sizeBeforeAdditions != refinementResult->addedGridPoints.size()) {
            std::cout << "Grid growth not correlated to refinement results (grid delta "
                      << sizeAfterAdditions - sizeBeforeAdditions << ", additions "
                      << refinementResult->addedGridPoints.size() << ")" << std::endl;
            exit(-1);
        }
    } else {
        // apply grid changes to the Cholesky factorization
        unsigned long currentGridVersion = learnerInstance->getLocalGridVersion(classIndex);
        if (!learnerInstance->checkGridStateConsistent(classIndex)) {
            std::cout << "Attempting to increment grid version on non consistent grid "
                      << classIndex
                      << " version " << currentGridVersion << std::endl;
            exit(-1);
        }

        learnerInstance->setLocalGridVersion(classIndex, currentGridVersion + 1);

        // Send class update in preparation for cholesky
        MPIMethods::sendRefinementUpdates(classIndex, refinementResult->deletedGridPointsIndexes,
                                          refinementResult->addedGridPoints);

        if (learnerInstance->getOffline()->isRefineable()) {
            AssignTaskResult assignTaskResult{};
            learnerInstance->getScheduler().assignTaskStaticTaskSize(
                    RECOMPUTE_CHOLESKY_DECOMPOSITION,
                    assignTaskResult);

            std::cout << "Assigning update of cholesky decomposition " << classIndex
                      << " to worker "
                      << assignTaskResult.workerID << std::endl;
            MPIMethods::assignSystemMatrixUpdate(assignTaskResult.workerID, classIndex);
        }
    }

    // update alpha vector
    D(unsigned long oldSize = densEst->getAlpha().size();)
    densEst->updateAlpha(&(refinementResult->deletedGridPointsIndexes),
                         refinementResult->addedGridPoints.size());
    D(std::cout << "Updated alpha vector " << classIndex << " (old size " << oldSize
                << ", new size " <<
                densEst->getAlpha().size() << ")" << std::endl;)
}

bool LearnerSGDEOnOffParallelHandler::checkRefinementNecessary(const std::string &refMonitor,
                                                               unsigned long refPeriod,
                                                               unsigned long totalInstances,
                                                               double currentValidError,
                                                               double currentTrainError,
                                                               unsigned long numberOfCompletedRefinements,
                                                               ConvergenceMonitor &monitor) {
    auto &offline = learnerInstance->getOffline();
    // access DBMatOnlineDE-objects of all classes in order
    // to apply adaptivity to the specific sparse grids later on

    // check if refinement should be performed
    if (refMonitor == "periodic") {
        // check periodic monitor
        if (offline->isRefineable() && (totalInstances > 0) && (totalInstances % refPeriod == 0) &&
            (numberOfCompletedRefinements < offline->getConfig().numRefinements_)) {
            return true;
        }
    } else if (refMonitor == "convergence") {
        // check convergence monitor
        if (learnerInstance->getValidationData() == nullptr) {
            throw sgpp::base::data_exception(
                    "No validation data for checking convergence provided!");
        }
        if (offline->isRefineable() &&
            (numberOfCompletedRefinements < offline->getConfig().numRefinements_)) {
            currentValidError = learnerInstance->getError(*learnerInstance->getValidationData());
            currentTrainError = learnerInstance->getError(
                    learnerInstance->getTrainData());  // if train dataset is large
            // use a subset for error
            // evaluation
            monitor.pushToBuffer(currentValidError, currentTrainError);
            if (monitor.nextRefCnt > 0) {
                monitor.nextRefCnt--;
            }
            if (monitor.nextRefCnt == 0) {
                return monitor.checkConvergence();
            }
        }
    }
    return false;
}

unsigned long
LearnerSGDEOnOffParallelHandler::handleDataAndZeroBasedRefinement(bool preCompute,
                                                                  MultiGridRefinementFunctor *func,
                                                                  unsigned long idx,
                                                                  base::Grid &grid,
                                                                  base::GridGenerator &gridGen) const {
    if (preCompute) {
        // precompute the evals (needs to be done once per step, before
        // any refinement is done
        func->preComputeEvaluations();
    }
    func->setGridIndex(idx);
    // perform refinement (zero-crossings-based / data-based)
    unsigned long gridSizeBeforeRefine = grid.getSize();
    gridGen.refine(*func);
    unsigned long gridSizeAfterRefine = grid.getSize();
    return gridSizeAfterRefine - gridSizeBeforeRefine;
}

unsigned long
LearnerSGDEOnOffParallelHandler::handleSurplusBasedRefinement(DBMatOnlineDE *densEst,
                                                              base::Grid &grid,
                                                              base::GridGenerator &gridGen) const {
    DataVector *alphaWork;  // required for surplus refinement
    // auxiliary variables
    DataVector p(learnerInstance->getTrainData().getDimension());


    std::unique_ptr<sgpp::base::OperationEval> opEval(op_factory::createOperationEval(grid));
    sgpp::base::HashGridStorage &gridStorage = grid.getStorage();
    alphaWork = &(densEst->getAlpha());
    DataVector alphaWeight(alphaWork->getSize());
    // determine surpluses
    for (unsigned long k = 0; k < gridStorage.getSize(); k++) {
        // sets values of p to the coordinates of the given GridPoint gp
        gridStorage.getPoint(k).getStandardCoordinates(p);
        // multiply k-th alpha with the evaluated function at grind-point
        // k
        alphaWeight[k] = alphaWork->get(k) * opEval->eval(*alphaWork, p);
    }

    // Perform Coarsening (surplus based)
    /*if (coarseCnt < maxCoarseNum) {
      HashCoarsening coarse_;
      //std::cout << "\n" << "Start coarsening\n";

      // Coarsening based on surpluses
      SurplusCoarseningFunctor scf(
        alphaWeight, coarseNumPoints, coarseThreshold);

      //std::cout << "Size before coarsening:" << grid->getSize() <<
    "\n";
      //int old_size = grid->getSize();
      coarse_.free_coarsen_NFirstOnly(
        grid->getStorage(), scf, alphaWeight, grid->getSize());

      std::cout << "Size after coarsening:" << grid->getSize() <<
    "\n\n";
      //int new_size = grid->getSize();

      deletedGridPoints.clear();
      deletedGridPoints = coarse_.getDeletedPoints();

      (*refineCoarse)[idx].first = deletedGridPoints;

      coarseCnt++;
    }*/

    // perform refinement (surplus based)
    unsigned long sizeBeforeRefine = grid.getSize();
    // simple refinement based on surpluses
    sgpp::base::SurplusRefinementFunctor srf(alphaWeight,
                                             learnerInstance->getOffline()->getConfig().ref_noPoints_);
    gridGen.refine(srf);
    unsigned long sizeAfterRefine = grid.getSize();
    return sizeAfterRefine - sizeBeforeRefine;
}

RefinementResult &LearnerSGDEOnOffParallelHandler::getRefinementResult(unsigned long classIndex) {
    return vectorRefinementResults[classIndex];
}

LearnerSGDEOnOffParallelHandler::LearnerSGDEOnOffParallelHandler(
        LearnerSGDEOnOffParallel *learnerInstance,
        size_t numClasses) {
    this->learnerInstance = learnerInstance;
    RefinementResult initResult{};
    vectorRefinementResults.insert(
            vectorRefinementResults.begin(), numClasses, initResult);

}
}  // namespace datadriven
}  // namespace sgpp


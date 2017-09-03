// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

//TODO: Delete this Definition after Refactoring
#define USE_GSL


#ifdef USE_GSL

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/algorithm/ConvergenceMonitor.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>
#include <sgpp/datadriven/application/work_in_progress/MPIMethods.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>


#include <chrono>
#include <thread>
#include <numeric>

using sgpp::base::Grid;
using sgpp::base::GridStorage;
using sgpp::base::GridGenerator;
using sgpp::base::OperationEval;
using sgpp::base::algorithm_exception;
using sgpp::base::data_exception;
using sgpp::base::SurplusRefinementFunctor;


namespace sgpp {
    namespace datadriven {

        LearnerSGDEOnOffParallel::LearnerSGDEOnOffParallel(sgpp::datadriven::DBMatDensityConfiguration &dconf,
                                                           Dataset &trainData, Dataset &testData,
                                                           Dataset *validationData,
                                                           sgpp::base::DataVector &classLabels, size_t numClassesInit,
                                                           bool usePrior, double beta, double lambda,
                                                           MPITaskScheduler &mpiTaskScheduler)
                : LearnerSGDEOnOff(dconf, trainData, testData, validationData, classLabels, numClassesInit, usePrior,
                                   beta, lambda), mpiTaskScheduler(mpiTaskScheduler) {

            RefinementResult initResult{};
            vectorRefinementResults.insert(vectorRefinementResults.begin(), getNumClasses(), initResult);
            localGridVersions.insert(localGridVersions.begin(), numClasses, 10);
            mpiTaskScheduler.setLearnerInstance(this);
            workerActive = true;

            MPIMethods::initMPI(this);
        }

        LearnerSGDEOnOffParallel::~LearnerSGDEOnOffParallel() {
            MPIMethods::finalizeMPI();
        }

        void LearnerSGDEOnOffParallel::train(size_t batchSize, size_t maxDataPasses,
                                             std::string refinementFunctorType, std::string refMonitor,
                                             size_t refPeriod, double accDeclineThreshold,
                                             size_t accDeclineBufferSize, size_t minRefInterval,
                                             bool enableCv, size_t nextCvStep) {

            if (!MPIMethods::isMaster()) {
                //TODO: Avoid queue size check
                while (workerActive || MPIMethods::getQueueSize() > 2) {
                    D(std::cout << "Client looping" << std::endl;)
                    MPIMethods::waitForAnyMPIRequestsToComplete();
                }
                std::cout << "Worker shutdown." << std::endl;
                MPIMethods::sendCommandNoArgs(MPI_MASTER_RANK, WORKER_SHUTDOWN_SUCCESS);
                //TODO: Avoid queue size check
                while (MPIMethods::getQueueSize() > 2) {
                    //TODO: Avoid queue size check
                    std::cout << "Waiting for " << MPIMethods::getQueueSize() - 2 << " operations to complete"
                              << std::endl;
                    MPIMethods::waitForAnyMPIRequestsToComplete();
                }
                std::cout << "Sent acknowledgement" << std::endl;
                return;
            }

            MPIMethods::processCompletedMPIRequests();

            // contains list of removed grid points and number of added grid points
            // (is updated in each refinement/coarsening step)
//            vectorRefinementResults = new ...

            // initialize counter for dataset passes
            size_t completedDataPasses = 0;

            // initialize refinement variables
            double currentValidError = 0.0;
            double currentTrainError = 0.0;
            // create convergence monitor object
            ConvergenceMonitor monitor{accDeclineThreshold, accDeclineBufferSize, minRefInterval};

            // counts number of performed refinements
            size_t numberOfCompletedRefinements = 0;

            // coarsening
            // size_t coarseCnt = 0;
            // size_t maxCoarseNum = 5;
            // size_t coarsePeriod = 50;
            // size_t coarseNumPoints = 1;
            // double coarseThreshold = 1.0;

            auto &onlineObjects = getDensityFunctions();

            size_t dim = getDimensionality();

            // pointer to the next batch (data points + class labels) to be processed
            Dataset dataBatch(batchSize, dim);

            // print initial grid size
            printGridSizeStatistics("#Initial grid size of grid ", onlineObjects);

            // auxiliary variable for accuracy (error) measurement
            // TODO Evil
//            double acc = getAccuracy();
//
//            avgErrors.append(1.0 - acc);

            // main loop which performs the training process
            while (completedDataPasses < maxDataPasses) {
                std::cout << "Start of data pass " << completedDataPasses << std::endl;

                std::cout << "#batch-size: " << batchSize << std::endl;

                // iterate over total number of batches
                while (processedPoints < trainData.getNumberInstances()) {
                    D(std::cout << "#processing batch: " << processedPoints << "\n";)

                    auto begin = std::chrono::high_resolution_clock::now();

                    // check if cross-validation should be performed
                    bool doCrossValidation = false;
                    if (enableCv) {
                        if (nextCvStep == processedPoints) {
                            doCrossValidation = true;
                            nextCvStep *= 5;
                        }
                    }

                    assignBatchToWorker(processedPoints, doCrossValidation);

                    processedPoints += dataBatch.getNumberInstances();

                    std::cout << processedPoints << " have already been assigned." << std::endl;


                    // Refinement only occurs on the Master Node

                    D(std::cout << "Checking if refinement is necessary." << std::endl;)
                    // check if refinement should be performed
                    if (checkRefinementNecessary(refMonitor, refPeriod, processedPoints, currentValidError,
                                                 currentTrainError, numberOfCompletedRefinements, monitor)) {

                        while (!checkReadyForRefinement()) {
                            D(std::cout << "Waiting for " << MPIMethods::getQueueSize()
                                        << " queue operations to complete before continuing" << std::endl;)
                            MPIMethods::waitForAnyMPIRequestsToComplete();
                        }

                        // if the Cholesky decomposition is chosen as factorization method
                        // refinement
                        // and coarsening methods can be applied

                        std::cout << "refinement at iteration: " << processedPoints << std::endl;
                        mpiTaskScheduler.onRefinementStarted();

                        doRefinementForAll(refinementFunctorType, refMonitor, onlineObjects, monitor);
                        numberOfCompletedRefinements += 1;
                        D(std::cout << "Refinement at " << processedPoints << " complete" << std::endl;)

                        //Send the grid component update
                        //Note: This was moved to updateClassVariablesAfterRefinement as it needs to run before the cholesky update
//                        MPIMethods::sendGridComponentsUpdate(vectorRefinementResults);
                    } else {
                        D(std::cout << "No refinement necessary" << std::endl;)
                    }

                    //TODO: Evil
//                    // save current error
//                    if ((processedPoints / batchSize) % 10 == 0) {
//                        acc = getAccuracy();
//                        std::cout << "Saving current error " << acc << std::endl;
//                        avgErrors.append(1.0 - acc);
//                    }

                    auto end = std::chrono::high_resolution_clock::now();
                    D(
                            std::cout << "Processing batch in "
                                      << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
                                      << "ms" << std::endl;
                            std::cout << "Processed " << processedPoints << " data points so far" << std::endl;
                    )
                }

                // Synchronize end of Data Pass
                std::cout << "End of data pass " << completedDataPasses << std::endl;

                completedDataPasses++;
                processedPoints = 0;
            }  // end while

            shutdown();

            std::cout << "#Training finished (This is MASTER)" << std::endl;

            //TODO This is also evil
//            error = 1.0 - getAccuracy();

        }

        bool LearnerSGDEOnOffParallel::checkReadyForRefinement() const {
            // All local grids in a consistent state and have OK from scheduler
            return mpiTaskScheduler.isReadyForRefinement()
                   && std::all_of(localGridVersions.begin(), localGridVersions.end(),
                                  [](size_t version) { return isVersionConsistent(version); });
        }

        size_t LearnerSGDEOnOffParallel::getDimensionality() {
            return trainData.getDimension();
        }

        void LearnerSGDEOnOffParallel::printGridSizeStatistics(const char *messageString,
                                                               ClassDensityContainer &onlineObjects) {
            // print initial grid size
            for (auto &onlineObject : onlineObjects) {
                auto densEst = onlineObject.first.get();
                Grid &grid = densEst->getOfflineObject().getGrid();
                std::cout << messageString << onlineObject.second << ", " << grid.getSize() << std::endl;

            }
        }

        void LearnerSGDEOnOffParallel::doRefinementForAll(const std::string &refinementFunctorType,
                                                          const std::string &refinementMonitorType,
                                                          const ClassDensityContainer &onlineObjects,
                                                          ConvergenceMonitor &monitor) {
            // acc = getAccuracy();
            // avgErrors.append(1.0 - acc);

            // bundle grids and surplus vector pointer needed for refinement
            // (for zero-crossings refinement, data-based refinement)
            std::vector<Grid *> grids;
            std::vector<DataVector *> alphas;
            for (size_t i = 0; i < getNumClasses(); i++) {
                auto densEst = onlineObjects[i].first.get();
                grids.push_back(&(densEst->getOfflineObject().getGrid()));
                alphas.push_back(&(densEst->getAlpha()));
            }
            bool levelPenalize = false;  // Multiplies penalzing term for fine levels
            bool preCompute = true;      // Precomputes and caches evals for zrcr
            MultiGridRefinementFunctor *func = nullptr;

            // Zero-crossing-based refinement
            ZeroCrossingRefinementFunctor funcZrcr{grids, alphas, offline->getConfig().ref_noPoints_,
                                                   levelPenalize, preCompute};

            // Data-based refinement. Needs a problem dependent coeffA. The values
            // can be determined by testing (aim at ~10 % of the training data is
            // to be marked relevant). Cross-validation or similar can/should be
            // employed
            // to determine this value.
            std::vector<double> coeffA;
            coeffA.push_back(1.2);  // ripley 1.2
            coeffA.push_back(1.2);  // ripley 1.2
            DataMatrix *trainDataRef = &(trainData.getData());
            DataVector *trainLabelsRef = &(trainData.getTargets());
            DataBasedRefinementFunctor funcData = DataBasedRefinementFunctor{
                    grids, alphas, trainDataRef, trainLabelsRef, offline->getConfig().ref_noPoints_,
                    levelPenalize, coeffA};

            if (refinementFunctorType == "zero") {
                func = &funcZrcr;
            } else if (refinementFunctorType == "data") {
                func = &funcData;
            }

            // perform refinement/coarsening for each grid
            for (size_t idx = 0; idx < this->getNumClasses(); idx++) {
                this->doRefinementForClass(refinementFunctorType,
                                           &(vectorRefinementResults[idx]),
                                           onlineObjects,
                                           preCompute, func, idx);
            }
            if (refinementMonitorType == "convergence") {
                monitor.nextRefCnt = monitor.minRefInterval;
            }

            //Wait for all new cholesky decompositions to come back
            MPIMethods::waitForIncomingMessageType(UPDATE_GRID, getNumClasses(), [](PendingMPIRequest &request) {
                auto *refinementResultNetworkMessage = (RefinementResultNetworkMessage *) request.buffer->payload;
                //Ensure it is a cholesky packet and the last in the sequence
                D(std::cout << "Test packet grid version " << refinementResultNetworkMessage->gridversion
                            << ", update type " << refinementResultNetworkMessage->updateType << std::endl;)
                return isVersionConsistent(refinementResultNetworkMessage->gridversion) &&
                       refinementResultNetworkMessage->updateType == CHOLESKY_DECOMPOSITION;
            });
        }

        void LearnerSGDEOnOffParallel::doRefinementForClass(const std::string &refType,
                                                            RefinementResult *refinementResult,
                                                            const ClassDensityContainer &onlineObjects,
                                                            bool preCompute,
                                                            MultiGridRefinementFunctor *refinementFunctor,
                                                            size_t classIndex) {
            // perform refinement/coarsening for grid which corresponds to current
            // index
            std::cout << "Refinement and coarsening for class: " << classIndex
                      << std::endl;
            auto densEst = onlineObjects[classIndex].first.get();
            Grid &grid = densEst->getOfflineObject().getGrid();

            size_t oldGridSize = grid.getSize();
            D(std::cout << "Size before adaptivity: " << oldGridSize
                        << std::endl;)

            base::GridGenerator &gridGen = grid.getGenerator();

            size_t numberOfNewPoints = 0;

            if (refType == "surplus") {
                numberOfNewPoints = handleSurplusBasedRefinement(densEst, grid, gridGen);

            } else if ((refType == "data") || (refType == "zero")) {
                numberOfNewPoints = handleDataAndZeroBasedRefinement(preCompute, refinementFunctor, classIndex, grid,
                                                                     gridGen);
            }

            size_t newGridSize = grid.getSize();
            std::cout << "grid size after adaptivity: " << newGridSize
                      << " (previously " << oldGridSize
                      << "), " << numberOfNewPoints << " new points on grid"
                      << std::endl;

            if (numberOfNewPoints != newGridSize - oldGridSize) {
                std::cout << "Reported grid sizes do not match up (refined " << numberOfNewPoints << ", old "
                          << oldGridSize << ", new " << newGridSize << ")" << std::endl;
                exit(-1);
            }

            D(std::cout << "Preparing refinement result update" << std::endl;)
            D(std::cout << "Clearing old added grid points" << std::endl;)
            refinementResult->addedGridPoints.clear();

            D(std::cout << "Clearing old deleted grid points" << std::endl;)
            refinementResult->deletedGridPointsIndexes.clear();


            size_t numDimensions = getDimensionality();
            //Collect new grid points into the refinement result for shipping
            for (size_t i = 0; i < numberOfNewPoints; i++) {
                LevelIndexVector levelIndexVector(numDimensions);
                base::HashGridPoint &gridPoint = grid.getStorage()[oldGridSize + i];
                for (size_t currentDimension = 0; currentDimension < numDimensions; currentDimension += 1) {
                    base::HashGridPoint::level_type pointLevel;
                    base::HashGridPoint::index_type pointIndex;
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

        void LearnerSGDEOnOffParallel::updateClassVariablesAfterRefinement(size_t classIndex,
                                                                           RefinementResult *refinementResult,
                                                                           DBMatOnlineDE *densEst) {

            base::Grid &grid = densEst->getOfflineObject().getGrid();

            if (!MPIMethods::isMaster()) {
                std::cout << "Applying refinement result class " << classIndex << " from master" << std::endl;
                D(std::cout << "Old grid " << classIndex << " size is " << grid.getSize() << std::endl;)


                size_t numDimensions = getDimensionality();

                //Delete the grid points removed on master thread
                grid.getStorage().deletePoints(refinementResult->deletedGridPointsIndexes);

                size_t sizeBeforeAdditions = grid.getSize();
                D(std::cout << "Grid size after deleting is " << sizeBeforeAdditions << std::endl;)

                //Add grid points added on master thread
                for (LevelIndexVector &levelIndexVector : refinementResult->addedGridPoints) {
                    size_t sizeBeforePoint = grid.getSize();
                    auto *gridPoint = new sgpp::base::HashGridStorage::point_type(numDimensions);
                    //TODO: What happens when other points are changed (ie Leaf boolean etc)

                    for (size_t currentDimension = 0; currentDimension < numDimensions; currentDimension++) {
                        gridPoint->set(currentDimension,
                                       levelIndexVector[currentDimension].level,
                                       levelIndexVector[currentDimension].index);
                    }
                    grid.getStorage().insert(*gridPoint);
                    D(
                            size_t sizeAfterPoint = grid.getSize();
                            if (sizeAfterPoint - sizeBeforePoint != 1) {
                                std::cout << "Inserted grid point but size change incorrect (old " << sizeBeforePoint
                                          << ", new " << sizeAfterPoint << "), point " << gridPoint->getHash()
                                          << std::endl;
                                printPoint(gridPoint);
                                exit(-1);
                            }
                    )
                }

                //TODO: This might be unnecessary
                grid.getStorage().recalcLeafProperty();

                size_t sizeAfterAdditions = grid.getSize();
                D(std::cout << "New grid " << classIndex << " size is " << sizeAfterAdditions << std::endl;)
                if (sizeAfterAdditions - sizeBeforeAdditions != refinementResult->addedGridPoints.size()) {
                    std::cout << "Grid growth not correlated to refinement results (grid delta "
                              << sizeAfterAdditions - sizeBeforeAdditions << ", additions "
                              << refinementResult->addedGridPoints.size() << ")" << std::endl;
                    exit(-1);
                }
            }
            else {
                // apply grid changes to the Cholesky factorization
                size_t currentGridVersion = getCurrentGridVersion(classIndex);
                if (!checkGridStateConsistent(classIndex)) {
                    std::cout << "Attempting to increment grid version on non consistent grid " << classIndex
                              << " version " << currentGridVersion << std::endl;
                    exit(-1);
                }

                setLocalGridVersion(classIndex, currentGridVersion + 1);

                // Send class update in preparation for cholesky
                MPIMethods::sendRefinementUpdates(classIndex, refinementResult->deletedGridPointsIndexes,
                                                  refinementResult->addedGridPoints);

                if (offline->isRefineable()) {
                    AssignTaskResult assignTaskResult{};
                    mpiTaskScheduler.assignTaskStaticTaskSize(RECOMPUTE_CHOLESKY_DECOMPOSITION, assignTaskResult);

                    std::cout << "Assigning update of cholesky decomposition " << classIndex << " to worker "
                              << assignTaskResult.workerID << std::endl;
                    MPIMethods::assignSystemMatrixUpdate(assignTaskResult.workerID, classIndex);
                }
            }

            // update alpha vector
            size_t oldSize = densEst->getAlpha().size();
            densEst->updateAlpha(&(refinementResult->deletedGridPointsIndexes),
                                 refinementResult->addedGridPoints.size());
            D(std::cout << "Updated alpha vector " << classIndex << " (old size " << oldSize << ", new size " <<
                        densEst->getAlpha().size() << ")" << std::endl;)
        }

        void
        LearnerSGDEOnOffParallel::computeNewCholeskyDecomposition(size_t classIndex, size_t gridversion) {

            // The first check is to ensure that all segments of an update have been received (intermediate segments set grid version to TEMPORARILY_INCONSISTENT)
            RefinementResult &refinementResult = vectorRefinementResults[classIndex];
            while (getCurrentGridVersion(classIndex) == GRID_TEMPORARILY_INCONSISTENT || (
                    refinementResult.deletedGridPointsIndexes.empty() &&
                    refinementResult.addedGridPoints.empty())) {
                D(std::cout << "Refinement results have not arrived yet (grid version "
                            << getCurrentGridVersion(classIndex)
                            << ", additions " << refinementResult.addedGridPoints.size() << ", deletions "
                            << refinementResult.deletedGridPointsIndexes.size() << "). Waiting..." << std::endl;)
                // Do not use waitForConsistent here, we want GRID_ADDITIONS or GRID_DELETIONS, not consistency

                MPIMethods::waitForIncomingMessageType(UPDATE_GRID);
                D(std::cout << "Updates have arrived. Attempting to resume." << std::endl;)
            }

            std::cout << "Computing cholesky modification for class " << classIndex << std::endl;

            DBMatOnlineDE *densEst = getDensityFunctions()[classIndex].first.get();
            auto &dbMatOfflineChol = dynamic_cast<DBMatOfflineChol &>(densEst->getOfflineObject());
            dbMatOfflineChol.choleskyModification(refinementResult.addedGridPoints.size(),
                                                  refinementResult.deletedGridPointsIndexes, densEst->getBestLambda());

            setLocalGridVersion(classIndex, gridversion);
            D(std::cout << "Send cholesky update to master for class " << classIndex << std::endl;)
            DataMatrix &newDecomposition = dbMatOfflineChol.getDecomposedMatrix();
            MPIMethods::sendCholeskyDecomposition(classIndex, newDecomposition, 0);
        }

        void LearnerSGDEOnOffParallel::assembleNextBatchData(Dataset *dataBatch, size_t *batchOffset) const {

            size_t batchSize = dataBatch->getNumberInstances();
            size_t dataDimensionality = dataBatch->getDimension();
            D(std::cout << "Assembling batch of size " << batchSize << " at offset " << *batchOffset << std::endl;)

            for (size_t j = 0; j < batchSize; j++) {
                base::DataVector dataPoint(dataDimensionality);
                trainData.getData().getRow(j + *batchOffset, dataPoint);
                double y = trainData.getTargets().get(j + *batchOffset);
                dataBatch->getData().setRow(j, dataPoint);
                dataBatch->getTargets().set(j, y);
            }
            *batchOffset += batchSize;

            D(std::cout << "Finished assembling batch" << std::endl;)
        }

        bool LearnerSGDEOnOffParallel::checkRefinementNecessary(const std::string &refMonitor, size_t refPeriod,
                                                                size_t totalInstances, double currentValidError,
                                                                double currentTrainError,
                                                                size_t numberOfCompletedRefinements,
                                                                ConvergenceMonitor &monitor) {

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
                if (validationData == nullptr) {
                    throw data_exception("No validation data for checking convergence provided!");
                }
                if (offline->isRefineable() && (numberOfCompletedRefinements < offline->getConfig().numRefinements_)) {
                    currentValidError = getError(*validationData);
                    currentTrainError = getError(trainData);  // if train dataset is large
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

        size_t
        LearnerSGDEOnOffParallel::handleDataAndZeroBasedRefinement(bool preCompute, MultiGridRefinementFunctor *func,
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

        size_t LearnerSGDEOnOffParallel::handleSurplusBasedRefinement(DBMatOnlineDE *densEst, base::Grid &grid,
                                                                      base::GridGenerator &gridGen) const {
            DataVector *alphaWork;  // required for surplus refinement
            // auxiliary variables
            DataVector p(trainData.getDimension());


            std::unique_ptr<OperationEval> opEval(op_factory::createOperationEval(grid));
            GridStorage &gridStorage = grid.getStorage();
            alphaWork = &(densEst->getAlpha());
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
            size_t sizeBeforeRefine = grid.getSize();
            // simple refinement based on surpluses
            SurplusRefinementFunctor srf(alphaWeight, offline->getConfig().ref_noPoints_);
            gridGen.refine(srf);
            size_t sizeAfterRefine = grid.getSize();
            return sizeAfterRefine - sizeBeforeRefine;
        }

        // Train from an entire Batch
        void LearnerSGDEOnOffParallel::train(Dataset &dataset, bool doCrossValidation) {
            size_t dim = dataset.getDimension();

            std::cout << "Starting train cycle (dataset size: " << dataset.getNumberInstances() << ")" << std::endl;

            // create an empty matrix for every class:
            std::vector<std::unique_ptr<DataMatrix>> classData;
            classData.reserve(getNumClasses());
            std::vector<std::pair<DataMatrix *, double>> trainDataClasses;
            trainDataClasses.reserve(getNumClasses());

            std::map<double, int> classIndices;  // maps class labels to indices

            D(std::cout << "Allocating class matrices" << std::endl;)
            allocateClassMatrices(dim, trainDataClasses, classIndices);

            D(std::cout << "Splitting batch into classes" << std::endl;)
            // split the data into the different classes:
            splitBatchIntoClasses(dataset, dim, trainDataClasses, classIndices);

            D(std::cout << "Computing density functions" << std::endl;)
            // compute density functions
            train(trainDataClasses, doCrossValidation);

            D(std::cout << "Finished train cycle." << std::endl;)
        }

        void
        LearnerSGDEOnOffParallel::splitBatchIntoClasses(const Dataset &dataset,
                                                        size_t dim,
                                                        const std::vector<std::pair<DataMatrix *, double>> &trainDataClasses,
                                                        std::map<double, int> &classIndices) const {
            // split the data into the different classes:
            for (size_t i = 0; i < dataset.getNumberInstances(); i++) {
                double classLabel = dataset.getTargets()[i];
                DataVector vec(dim);
                dataset.getData().getRow(i, vec);
                auto &classPairDataMatrixDouble = trainDataClasses[classIndices[classLabel]];
                classPairDataMatrixDouble.first->appendRow(vec);
            }
        }

        void LearnerSGDEOnOffParallel::allocateClassMatrices(size_t dim,
                                                             std::vector<std::pair<base::DataMatrix *, double>> &trainDataClasses,
                                                             std::map<double, int> &classIndices) const {
            int index = 0;
            for (size_t i = 0; i < getNumClasses(); i++) {
                auto *m = new base::DataMatrix(0, dim);
                std::pair<base::DataMatrix *, double> p(m, classLabels[i]);
                trainDataClasses.push_back(p);
                classIndices.insert(std::pair<double, int>(classLabels[i], index));
                index++;
            }
        }

        //Train from a Batch already split up into its classes
        void
        LearnerSGDEOnOffParallel::train(std::vector<std::pair<sgpp::base::DataMatrix *, double> > &trainDataClasses,
                                        bool doCrossValidation) {

            //Calculate the total number of data points
            size_t numberOfDataPoints = 0;
            for (auto &trainDataClass : trainDataClasses) {
                numberOfDataPoints += trainDataClass.first->getSize();
            }

            // Learn from each Class
            for (size_t i = 0; i < trainDataClasses.size(); i++) {
                std::pair<sgpp::base::DataMatrix *, double> p = trainDataClasses[i];

                if ((*p.first).getNrows() > 0) {
                    // update density function for current class
                    RefinementResult &classRefinementResult = vectorRefinementResults[i];
                    std::cout << "Calling compute density function class " << i << " (refinement +"
                              << classRefinementResult.addedGridPoints.size() << ", -"
                              << classRefinementResult.deletedGridPointsIndexes.size() << ")" << std::endl;
                    densityFunctions[i].first->computeDensityFunction(
                            *p.first, true, doCrossValidation, &classRefinementResult.deletedGridPointsIndexes,
                            classRefinementResult.addedGridPoints.size());
                    D(std::cout << "Clearing the refinement results class " << i << std::endl;)
                    classRefinementResult.deletedGridPointsIndexes.clear();
                    classRefinementResult.addedGridPoints.clear();

                    if (usePrior) {
                        //TODO: This probably doesn't work
                        double newPrior = ((this->prior[p.second] * static_cast<double>(processedPoints)) +
                                           static_cast<double>(p.first->getSize())) /
                                          (static_cast<double>(numberOfDataPoints) +
                                           static_cast<double>(processedPoints));
                        D(std::cout << "Setting prior[" << p.second << "] to " << newPrior << std::endl;)
                        this->prior[p.second] =
                                newPrior;
                    } else {
                        D(std::cout << "Setting prior[" << p.second << "] to 1.0" << std::endl;)
                        this->prior[p.second] = 1.;
                    }
                }
            }

            this->processedPoints += numberOfDataPoints;
            trained = true;
        }

        void LearnerSGDEOnOffParallel::shutdown() {
            if (MPIMethods::isMaster()) {
                std::cout << "Broadcasting shutdown" << std::endl;
                MPIMethods::bcastCommandNoArgs(SHUTDOWN);
                MPIMethods::waitForIncomingMessageType(WORKER_SHUTDOWN_SUCCESS, MPIMethods::getWorldSize() - 1);
            } else {
                workerActive = false;
            }
        }

        void LearnerSGDEOnOffParallel::workBatch(Dataset dataset, size_t batchOffset, bool doCrossValidation) {


            waitForAllGridsConsistent();

            // assemble next batch
            std::cout << "Learning with batch of size " << dataset.getNumberInstances()
                      << " at offset " << batchOffset << std::endl;
            assembleNextBatchData(&dataset, &batchOffset);
            D(std::cout << "Batch of size " << dataset.getNumberInstances() << " assembled, starting with training."
                        << std::endl;)

            // train the model with current batch
            train(dataset, doCrossValidation);

            // Batch offset was already modified by assembleNextBatch
            D(std::cout << "Batch " << batchOffset - dataset.getNumberInstances() << " completed." << std::endl;)
            auto &densityFunctions = getDensityFunctions();
            for (size_t classIndex = 0; classIndex < getNumClasses(); classIndex++) {
                D(std::cout << "Sending alpha values to master for class " << classIndex << " with grid version "
                            << getCurrentGridVersion(classIndex) << std::endl;)
                auto &classDensityContainer = densityFunctions[classIndex];
                DataVector alphaVector = classDensityContainer.first->getAlpha();
                MPIMethods::sendMergeGridNetworkMessage(classIndex, batchOffset, dataset.getNumberInstances(),
                                                        alphaVector);

                DataVector &dataVector = getDensityFunctions()[classIndex].first->getAlpha();
                D(std::cout << "Local alpha sum " << classIndex << " is now "
                            << std::accumulate(dataVector.begin(), dataVector.end(), 0.0) << std::endl;)

            }
            D(std::cout << "Completed work batch " << batchOffset - dataset.getNumberInstances()
                        << " requested by master." << std::endl;)
        }

        void LearnerSGDEOnOffParallel::waitForAllGridsConsistent() {
            size_t classIndex = 0;
            while (classIndex < localGridVersions.size()) {
                if (!checkGridStateConsistent(classIndex)) {
                    std::cout << "Attempted to train from an inconsistent grid "
                              << classIndex << " version " << getCurrentGridVersion(classIndex) << std::endl;
                    MPIMethods::waitForGridConsistent(classIndex);
                    //start over, waiting might have changed other grids
                    classIndex = 0;
                } else {
                    classIndex++;
                }
            }
        }

        bool LearnerSGDEOnOffParallel::isVersionConsistent(size_t version) { return version >= 10; }

        size_t
        LearnerSGDEOnOffParallel::assignBatchToWorker(size_t batchOffset, bool doCrossValidation) {
            AssignTaskResult assignTaskResult{};
            mpiTaskScheduler.assignTaskVariableTaskSize(TRAIN_FROM_BATCH, assignTaskResult);

            if (assignTaskResult.taskSize + batchOffset > trainData.getNumberInstances()) {
                std::cout << "Shortening last batch." << std::endl;
                assignTaskResult.taskSize = trainData.getNumberInstances() - batchOffset;
            }

            std::cout << "Assigning batch " << batchOffset
                      << " to worker " << assignTaskResult.workerID
                      << " with size " << assignTaskResult.taskSize << std::endl;
            MPIMethods::assignBatch(assignTaskResult.workerID, batchOffset, assignTaskResult.taskSize,
                                    doCrossValidation);
            return assignTaskResult.taskSize;
        }


        void LearnerSGDEOnOffParallel::mergeAlphaValues(unsigned long classIndex, unsigned long gridVersion,
                                                        DataVector dataVector, unsigned long batchOffset,
                                                        unsigned long batchSize, bool isLastPacketInSeries) {
            MPIMethods::waitForGridConsistent(classIndex);


            D(
                    std::cout << "Remote alpha sum " << classIndex << " is "
                        << std::accumulate(dataVector.begin(), dataVector.end(), 0.0) << std::endl;
                    std::cout << "Batch size is " << batchSize << std::endl;
            )

            if (!isVersionConsistent(gridVersion)) {
                std::cout << "Received merge request with inconsistent grid " << classIndex << " version "
                          << gridVersion << std::endl;
                exit(-1);
            }

            size_t localGridVersion = getCurrentGridVersion(classIndex);
            if (isLastPacketInSeries) {
                mpiTaskScheduler.onMergeRequestIncoming(batchOffset, batchSize, gridVersion, localGridVersion);
            }

            if (gridVersion != localGridVersion) {
                D(std::cout << "Received merge grid request with incorrect grid version!"
                            << " local: " << localGridVersion
                            << ", remote: " << gridVersion
                            << std::endl;)
                if (gridVersion + 1 == localGridVersion) {
                    RefinementResult &refinementResult = vectorRefinementResults[classIndex];
                    std::list<size_t> &deletedPoints = refinementResult.deletedGridPointsIndexes;
                    std::list<LevelIndexVector> &addedPoints = refinementResult.addedGridPoints;

                    D(std::cout << "Attempting to automatically compensate for outdated grid." << std::endl <<
                                "Refinement result has " << addedPoints.size() << " additions and "
                                << deletedPoints.size() << " deletions" << std::endl <<
                                "The original remote size is " << dataVector.size() << std::endl;)
                    if (!deletedPoints.empty()
                        || !addedPoints.empty()) {
                        D(std::cout << "Found necessary refinement data" << std::endl;)

                        //See DBMatOnlineDe::updateAlpha()
                        if (!deletedPoints.empty()) {
                            D(std::cout << "Copying vector (deleting deleted grid points)." << std::endl;)
                            DataVector newAlpha{dataVector.getSize() - deletedPoints.size() + addedPoints.size()};
                            for (size_t i = 0; i < dataVector.getSize(); i++) {
                                if (std::find(deletedPoints.begin(), deletedPoints.end(), i) != deletedPoints.end()) {
                                    continue;
                                }

                                newAlpha.append(dataVector.get(i));
                            }
                            // set new alpha
                            dataVector = newAlpha;
                        }
                        dataVector.resizeZero(dataVector.size() + addedPoints.size());
                        D(std::cout << "New alpha vector is now " << dataVector.size() << " elements long."
                                    << std::endl;)
                    } else {
                        std::cout << "Refinement data has already been deleted." << std::endl
                                  << "This is probably because the master is training from batches. " << std::endl
                                  << "Cannot compensate, will now fail." << std::endl;
                        exit(-1);
                    }
                } else {
                    std::cout << "Merge request " << batchOffset << ", size " << batchSize
                              << ", older than one refinement cycle. Increase the refinement period."
                              << std::endl;
                    exit(-1);
                }

            }


            DataVector &localAlpha = getDensityFunctions()[classIndex].first->getAlpha();
            if (localAlpha.size() != dataVector.size()) {
                std::cout << "Received merge request with incorrect size (local " << localAlpha.size() << ", remote "
                          << dataVector.size() << "), local version is " << localGridVersions[classIndex]
                          << std::endl;
                exit(-1);
            }

            if (usePrior) {
                //TODO: Not implemented
                throw algorithm_exception("Use prior not implemented");
            } else {
                D(std::cout << "Setting prior [" << classLabels[classIndex] << "] to -1" << std::endl;)
                prior[classLabels[classIndex]] = 1.0;
            }

            D(std::cout << "Local alpha sum " << classIndex << " was "
                        << std::accumulate(dataVector.begin(), dataVector.end(), 0.0) << std::endl;)
            localAlpha.add(dataVector);
            D(std::cout << "Local alpha sum " << classIndex << " is now "
                        << std::accumulate(dataVector.begin(), dataVector.end(), 0.0) << std::endl;)
        }

        size_t LearnerSGDEOnOffParallel::getCurrentGridVersion(size_t classIndex) {
            return localGridVersions[classIndex];
        }

        bool LearnerSGDEOnOffParallel::checkGridStateConsistent(size_t classIndex) {
            if (localGridVersions.size() < classIndex) {
                throw algorithm_exception("Received request for consistency of class " + classIndex);
            }
            return isVersionConsistent(localGridVersions[classIndex]);
        }

        void LearnerSGDEOnOffParallel::setLocalGridVersion(size_t classIndex, size_t gridVersion) {
            D(
                    if (localGridVersions[classIndex] != gridVersion) {
                        std::cout << "Grid " << classIndex << " now has version " << gridVersion << " (previously "
                                  << localGridVersions[classIndex] << ")" << std::endl;
                    }
            )
            if (!checkGridStateConsistent(classIndex) && isVersionConsistent(gridVersion)) {
                std::cout << "Grid " << classIndex << " has been fully updated to version " << gridVersion << std::endl;
            }
            localGridVersions[classIndex] = gridVersion;
        }

        RefinementResult &LearnerSGDEOnOffParallel::getRefinementResult(size_t classIndex) {
            return vectorRefinementResults[classIndex];
        }

        void LearnerSGDEOnOffParallel::printPoint(base::HashGridStorage::point_type *gridPoint) {
            std::cout << "Grid point " << gridPoint->getHash() << std::endl;

            for (size_t currentDimension = 0; currentDimension < getDimensionality(); currentDimension++) {
                sgpp::base::HashGridPoint::level_type level = 0;
                sgpp::base::HashGridPoint::index_type index = 0;
                gridPoint->get(currentDimension, level, index);
                std::cout << "Dimension " << currentDimension
                          << ", level " << level << ", index "
                          << index << std::endl;
            }

        }

    }  // namespace datadriven
}  // namespace sgpp

#endif /* USE_GSL */

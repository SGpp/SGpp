// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/learnersgdeonoffparallel/LearnerSGDEOnOffParallel.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp>
#include <sgpp/datadriven/algorithm/GridFactory.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorConvergence.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/MPIMethods.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/MPITaskScheduler.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/RefinementHandler.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <climits>
#include <list>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using sgpp::base::algorithm_exception;
using sgpp::base::data_exception;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridStorage;
using sgpp::base::OperationEval;
using sgpp::base::SurplusRefinementFunctor;

namespace sgpp {
namespace datadriven {

static const int MINIMUM_CONSISTENT_GRID_VERSION = 10;

LearnerSGDEOnOffParallel::LearnerSGDEOnOffParallel(
    sgpp::base::RegularGridConfiguration &gridConfig,
    sgpp::base::AdaptivityConfiguration &adaptivityConfig,
    sgpp::datadriven::RegularizationConfiguration &regularizationConfig,
    sgpp::datadriven::DensityEstimationConfiguration &densityEstimationConfig, Dataset &trainData,
    Dataset &testData, Dataset *validationData, sgpp::base::DataVector &classLabels,
    size_t numClassesInit, bool usePrior, double beta, MPITaskScheduler &mpiTaskScheduler)
    : trainData{trainData},
      testData{testData},
      validationData{validationData},
      classLabels{classLabels},
      numClasses{numClassesInit},
      usePrior{usePrior},
      prior{},
      beta{beta},
      trained{false},
      offlineContainer{},
      densityFunctions{},
      processedPoints{0},
      avgErrors{0},
      mpiTaskScheduler(mpiTaskScheduler),
      refinementHandler(nullptr, 0),
      gridConfig{gridConfig},
      adaptivityConfig{adaptivityConfig},
      regularizationConfig{regularizationConfig},
      densityEstimationConfig{densityEstimationConfig} {
  // initialize offline object
  // Create an offline object that serves as template for all classes
  offline = std::unique_ptr<DBMatOffline>{DBMatOfflineFactory::buildOfflineObject(
      gridConfig, adaptivityConfig, regularizationConfig, densityEstimationConfig)};

  // Create a grid TODO(fuchsgruber): Maybe remove the learner class entirely??
  grids.reserve(numClasses);
  densityFunctions.reserve(numClasses);
  offlineContainer.reserve(numClasses);
  alphas.reserve(numClasses);
  GridFactory gridFactory;
  for (size_t classIndex = 0; classIndex < numClasses; classIndex++) {
    // Create a grid
    std::unique_ptr<Grid> grid =
        std::unique_ptr<Grid>{gridFactory.createGrid(gridConfig, std::set<std::set<size_t>>())};
    std::unique_ptr<DBMatOffline> offlineCloned = std::unique_ptr<DBMatOffline>{offline->clone()};
    offlineCloned->buildMatrix(grid.get(), regularizationConfig);
    offlineCloned->decomposeMatrix(regularizationConfig, densityEstimationConfig);
    auto densEst = std::unique_ptr<DBMatOnlineDE>{DBMatOnlineDEFactory::buildDBMatOnlineDE(
        *offlineCloned, *grid, regularizationConfig.lambda_, beta)};
    densityFunctions.emplace_back(std::make_pair(std::move(densEst), classIndex));
    prior.emplace(classLabels[classIndex], 0.0);
    DataVector *alpha = new DataVector(offlineCloned->getGridSize());
    alphas.emplace_back(alpha);
    grids.emplace_back(std::move(grid));
    offlineContainer.emplace_back(std::move(offlineCloned));
  }

  localGridVersions.insert(localGridVersions.begin(), numClasses, MINIMUM_CONSISTENT_GRID_VERSION);
  mpiTaskScheduler.setLearnerInstance(this);
  workerActive = true;

  refinementHandler = RefinementHandler(this, numClassesInit);

  MPIMethods::initMPI(this);
}

LearnerSGDEOnOffParallel::~LearnerSGDEOnOffParallel() { MPIMethods::finalizeMPI(); }

Grid &LearnerSGDEOnOffParallel::getGrid(size_t classIndex) { return *(grids[classIndex]); }

size_t LearnerSGDEOnOffParallel::getNumClasses() const { return numClasses; }

double LearnerSGDEOnOffParallel::getAccuracy() const {
  DataVector computedLabels(testData.getNumberInstances());
  predict(testData.getData(), computedLabels);
  size_t correct = 0;
  size_t correctLabel1 = 0;
  size_t correctLabel2 = 0;
  for (size_t i = 0; i < computedLabels.getSize(); i++) {
    if (computedLabels.get(i) == testData.getTargets().get(i)) {
      correct++;
    }
    if ((computedLabels.get(i) == -1.0) && (testData.getTargets().get(i) == -1.0)) {
      correctLabel1++;
    } else if ((computedLabels.get(i) == 1.0) && (testData.getTargets().get(i) == 1.0)) {
      correctLabel2++;
    }
  }
  // std::cout << "correct label (-1): " << correctLabel1 << "\n";
  // std::cout << "correct label (1): " << correctLabel2 << "\n";

  double acc = static_cast<double>(correct) / static_cast<double>(computedLabels.getSize());
  return acc;
}

void LearnerSGDEOnOffParallel::predict(DataMatrix &data, DataVector &result) const {
  // calculate per class densities
  std::vector<DataVector> perClassDensities;
  for (auto &densityFunction : densityFunctions) {
    perClassDensities.emplace_back(data.getNrows());
    size_t classIndex = densityFunction.second;
    densityFunction.first->eval(*(alphas[classIndex]), data, perClassDensities.back(),
                                *(grids[classIndex]), true);
    perClassDensities.back().mult(prior.at(classLabels[classIndex]));
  }

// now select the appropriate class
#pragma omp parallel for
  for (auto point = 0u; point < data.getNrows(); point++) {
    auto bestClass = 0.0;
    auto maxDensity = std::numeric_limits<double>::max() * (-1);
    for (auto classNum = 0u; classNum < numClasses; classNum++) {
      auto density = perClassDensities[classNum][point];
      if (density > maxDensity) {
        maxDensity = density;
        bestClass = classLabels[densityFunctions[classNum].second];
      }
    }
    if (bestClass == 0) {
      std::cerr << "LearnerSGDEOnOff: Warning: no best class found!\n";
    }
    result[point] = bestClass;
  }
}

double LearnerSGDEOnOffParallel::getError(Dataset &dataset) const {
  double res = -1.0;

  DataVector computedLabels(dataset.getNumberInstances());
  predict(dataset.getData(), computedLabels);
  size_t correct = 0;
  for (size_t i = 0; i < computedLabels.getSize(); i++) {
    if (computedLabels.get(i) == dataset.getTargets().get(i)) {
      correct++;
    }
  }

  double acc = static_cast<double>(correct) / static_cast<double>(computedLabels.getSize());
  res = 1.0 - acc;

  return res;
}

void LearnerSGDEOnOffParallel::updateAlpha(size_t classIndex, std::list<size_t> *deletedPoints,
                                           size_t newPoints) {
  DataVector &alpha = *(alphas[classIndex]);
  if (alpha.getSize() != 0 && deletedPoints != nullptr && !deletedPoints->empty()) {
    std::vector<size_t> vecDeletedPoints{std::begin(*deletedPoints), std::end(*deletedPoints)};
    alpha.remove(vecDeletedPoints);
  }
  if (newPoints > 0) {
    alpha.resizeZero(alpha.getSize() + newPoints);
  }
}

void LearnerSGDEOnOffParallel::trainParallel(size_t batchSize, size_t maxDataPasses,
                                             std::string refinementFunctorType,
                                             std::string refMonitor, size_t refPeriod,
                                             double accDeclineThreshold,
                                             size_t accDeclineBufferSize, size_t minRefInterval) {
  if (!MPIMethods::isMaster()) {
    while (workerActive || MPIMethods::hasPendingOutgoingRequests()) {
      D(std::cout << "Client looping" << std::endl;)
      MPIMethods::waitForAnyMPIRequestsToComplete();
    }
    std::cout << "Worker shutdown." << std::endl;
    MPIMethods::sendCommandNoArgs(MPI_MASTER_RANK, WORKER_SHUTDOWN_SUCCESS);
    while (MPIMethods::hasPendingOutgoingRequests()) {
      std::cout << "Waiting for all outgoing requests to complete" << std::endl;
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
  RefinementMonitorConvergence monitor{accDeclineThreshold, accDeclineBufferSize, minRefInterval};

  // counts number of performed refinements
  size_t numberOfCompletedRefinements = 0;

  // coarsening
  // size_t coarseCnt = 0;
  // size_t maxCoarseNum = 5;
  // size_t coarsePeriod = 50;
  // size_t coarseNumPoints = 1;
  // double coarseThreshold = 1.0;

  auto &onlineObjects = getDensityFunctions();

  // print initial grid size
  printGridSizeStatistics("#Initial grid size of grid ", onlineObjects);

  // main loop which performs the training process
  while (completedDataPasses < maxDataPasses) {
    std::cout << "Start of data pass " << completedDataPasses << std::endl;

    std::cout << "#batch-size: " << batchSize << std::endl;

    // iterate over total number of batches
    while (processedPoints < trainData.getNumberInstances()) {
      D(std::cout << "#processing batch: " << processedPoints << "\n";
        auto begin = std::chrono::high_resolution_clock::now();)

      size_t batchSize = assignBatchToWorker(processedPoints, false);
      processedPoints += batchSize;

      std::cout << processedPoints << " have already been assigned." << std::endl;

      // Refinement only occurs on the Master Node

      D(std::cout << "Checking if refinement is necessary." << std::endl;)
      // check if refinement should be performed
      size_t refinementsNecessary = refinementHandler.checkRefinementNecessary(
          refMonitor, refPeriod, batchSize, currentValidError, currentTrainError,
          numberOfCompletedRefinements, monitor, adaptivityConfig);
      while (refinementsNecessary > 0) {
        while (!refinementHandler.checkReadyForRefinement()) {
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
        refinementsNecessary--;
        D(std::cout << "Refinement at " << processedPoints << " complete" << std::endl;)

        // Send the grid component update
        // Note: This was moved to updateClassVariablesAfterRefinement
        // as it needs to run before the system matrix update
        //        MPIMethods::sendGridComponentsUpdate(vectorRefinementResults);
      }

      D(auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Processing batch in "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
                  << "ms" << std::endl;
        std::cout << "Processed " << processedPoints << " data points so far" << std::endl;)
    }

    // Synchronize end of Data Pass
    std::cout << "End of data pass " << completedDataPasses << std::endl;

    completedDataPasses++;
    processedPoints = 0;
  }  // end while

  shutdownMPINodes();

  std::cout << "#Training finished (This is MASTER)" << std::endl;
}

size_t LearnerSGDEOnOffParallel::getDimensionality() { return trainData.getDimension(); }

std::vector<std::pair<std::unique_ptr<DBMatOnlineDE>, size_t>>
    &LearnerSGDEOnOffParallel::getDensityFunctions() {
  return densityFunctions;
}

void LearnerSGDEOnOffParallel::printGridSizeStatistics(
    const char *messageString,
    std::vector<std::pair<std::unique_ptr<DBMatOnlineDE>, size_t>> &onlineObjects) {
  // print initial grid size
  for (size_t idx = 0; idx < onlineObjects.size(); idx++) {
    auto &onlineObject = onlineObjects[idx];
    Grid &grid = *(grids[idx]);
    std::cout << messageString << onlineObject.second << ", " << grid.getSize() << std::endl;
  }
}

void LearnerSGDEOnOffParallel::doRefinementForAll(
    const std::string &refinementFunctorType, const std::string &refinementMonitorType,
    const std::vector<std::pair<std::unique_ptr<DBMatOnlineDE>, size_t>> &onlineObjects,
    RefinementMonitor &monitor) {
  // acc = getAccuracy();
  // avgErrors.append(1.0 - acc);

  // bundle grids and surplus vector pointer needed for refinement
  // (for zero-crossings refinement, data-based refinement)
  bool levelPenalize = false;  // Multiplies penalzing term for fine levels
  bool preCompute = true;      // Precomputes and caches evals for zrcr
  MultiGridRefinementFunctor *func = nullptr;

  // Copy the vector of grids for typecasting
  std::vector<Grid *> gridVector;
  std::vector<double> priorVector;
  gridVector.reserve(grids.size());
  for (size_t idx = 0; idx < grids.size(); idx++) {
    gridVector.push_back(&(*(grids[idx])));
    priorVector.push_back(prior.at(classLabels[idx]));
  }

  // Zero-crossing-based refinement
  ZeroCrossingRefinementFunctor funcZrcr{gridVector,    alphas,
                                         priorVector,   adaptivityConfig.numRefinementPoints_,
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
      gridVector,    alphas,         priorVector,
      trainDataRef,  trainLabelsRef, adaptivityConfig.numRefinementPoints_,
      levelPenalize, coeffA};

  if (refinementFunctorType == "zero") {
    func = &funcZrcr;
  } else if (refinementFunctorType == "data") {
    func = &funcData;
  }

  // perform refinement/coarsening for each grid
  for (size_t classIndex = 0; classIndex < this->getNumClasses(); classIndex++) {
    refinementHandler.doRefinementForClass(
        refinementFunctorType, &(refinementHandler.getRefinementResult(classIndex)), onlineObjects,
        *(grids[classIndex]), *(alphas[classIndex]), preCompute, func, classIndex,
        adaptivityConfig);
  }

  // Wait for all new system matrix decompositions to come back
  CHECK_SIZE_T_TO_INT(getNumClasses())
  MPIMethods::waitForIncomingMessageType(
      UPDATE_GRID, static_cast<int>(getNumClasses()), [](PendingMPIRequest &request) {
        auto *refinementResultNetworkMessage = static_cast<RefinementResultNetworkMessage *>(
            static_cast<void *>(request.buffer->payload));
        // Ensure it is a system matrix packet and the last in the sequence
        D(std::cout << "Test packet grid version " << refinementResultNetworkMessage->gridversion
                    << ", update type " << refinementResultNetworkMessage->updateType << std::endl;)
        return isVersionConsistent(refinementResultNetworkMessage->gridversion) &&
               refinementResultNetworkMessage->updateType == SYSTEM_MATRIX_DECOMPOSITION;
      });
}

void LearnerSGDEOnOffParallel::computeNewSystemMatrixDecomposition(size_t classIndex,
                                                                   size_t gridVersion) {
  // The first check is to ensure that all segments of an update have been
  // received
  // (intermediate segments set grid version to TEMPORARILY_INCONSISTENT)
  RefinementResult &refinementResult = refinementHandler.getRefinementResult(classIndex);
  while (getLocalGridVersion(classIndex) != GRID_RECEIVED_ADDED_POINTS ||
         (refinementResult.deletedGridPointsIndices.empty() &&
          refinementResult.addedGridPoints.empty())) {
    D(std::cout << "Refinement results have not arrived yet (grid version "
                << getLocalGridVersion(classIndex) << ", additions "
                << refinementResult.addedGridPoints.size() << ", deletions "
                << refinementResult.deletedGridPointsIndices.size() << "). Waiting..."
                << std::endl;)
    // Do not use waitForConsistent here, we want GRID_ADDITIONS or
    // GRID_DELETIONS, not consistency

    MPIMethods::waitForIncomingMessageType(UPDATE_GRID);
    D(std::cout << "Updates have arrived. Attempting to resume." << std::endl;)
  }

  std::cout << "Computing system matrix modification for class " << classIndex << "(+"
            << refinementResult.addedGridPoints.size() << ", -"
            << refinementResult.deletedGridPointsIndices.size() << ")" << std::endl;

  DBMatOnlineDE *densEst = getDensityFunctions()[classIndex].first.get();
  DBMatOffline &dbMatOffline = densEst->getOfflineObject();

  std::vector<size_t> idxToDelete{std::begin(refinementResult.deletedGridPointsIndices),
                                  std::end(refinementResult.deletedGridPointsIndices)};
  densEst->updateSystemMatrixDecomposition(densityEstimationConfig, *(grids[classIndex]),
                                           refinementResult.addedGridPoints.size(), idxToDelete,
                                           regularizationConfig.lambda_);

  setLocalGridVersion(classIndex, gridVersion);
  D(std::cout << "Send system matrix update to master for class " << classIndex << std::endl;)
  DataMatrix &newDecomposition = dbMatOffline.getDecomposedMatrix();
  MPIMethods::sendSystemMatrixDecomposition(classIndex, newDecomposition, 0);
}

void LearnerSGDEOnOffParallel::assembleNextBatchData(Dataset *dataBatch,
                                                     size_t *batchOffset) const {
  size_t batchSize = dataBatch->getNumberInstances();
  size_t dataDimensionality = dataBatch->getDimension();
  D(std::cout << "Assembling batch of size " << batchSize << " at offset " << *batchOffset
              << std::endl;)

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

// Train from an entire Batch
void LearnerSGDEOnOffParallel::train(Dataset &dataset, bool doCrossValidation) {
  size_t dim = dataset.getDimension();

  std::cout << "Starting train cycle (dataset size: " << dataset.getNumberInstances() << ")"
            << std::endl;

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

void LearnerSGDEOnOffParallel::splitBatchIntoClasses(
    const Dataset &dataset, size_t dim,
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

void LearnerSGDEOnOffParallel::allocateClassMatrices(
    size_t dim, std::vector<std::pair<DataMatrix *, double>> &trainDataClasses,
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

// Train from a Batch already split up into its classes
void LearnerSGDEOnOffParallel::train(std::vector<std::pair<DataMatrix *, double>> &trainDataClasses,
                                     bool doCrossValidation) {
  // Calculate the total number of data points
  size_t numberOfDataPoints = 0;
  for (auto &trainDataClass : trainDataClasses) {
    numberOfDataPoints += trainDataClass.first->getSize();
  }

  // Learn from each Class
  for (size_t classIndex = 0; classIndex < trainDataClasses.size(); classIndex++) {
    std::pair<sgpp::base::DataMatrix *, double> p = trainDataClasses[classIndex];

    if ((*p.first).getNrows() > 0) {
      // update density function for current class
      RefinementResult &classRefinementResult = refinementHandler.getRefinementResult(classIndex);
      std::cout << "Calling compute density function class " << classIndex << " (refinement +"
                << classRefinementResult.addedGridPoints.size() << ", -"
                << classRefinementResult.deletedGridPointsIndices.size() << ")" << std::endl;
      densityFunctions[classIndex].first->computeDensityFunction(
          *(alphas[classIndex]), *p.first, *(grids[classIndex]), densityEstimationConfig, true,
          doCrossValidation);
      D(std::cout << "Clearing the refinement results class " << classIndex << std::endl;)
      classRefinementResult.deletedGridPointsIndices.clear();
      classRefinementResult.addedGridPoints.clear();

      if (usePrior) {
        double newPrior =
            ((this->prior[p.second] * static_cast<double>(processedPoints)) +
             static_cast<double>(p.first->getSize())) /
            (static_cast<double>(numberOfDataPoints) + static_cast<double>(processedPoints));
        D(std::cout << "Setting prior[" << p.second << "] to " << newPrior << std::endl;)
        this->prior[p.second] = newPrior;
      } else {
        D(std::cout << "Setting prior[" << p.second << "] to 1.0" << std::endl;)
        this->prior[p.second] = 1.;
      }
    }
  }

  this->processedPoints += numberOfDataPoints;
  trained = true;
}

void LearnerSGDEOnOffParallel::shutdownMPINodes() {
  if (MPIMethods::isMaster()) {
    std::cout << "Broadcasting shutdown" << std::endl;
    MPIMethods::bcastCommandNoArgs(SHUTDOWN);
    MPIMethods::waitForIncomingMessageType(WORKER_SHUTDOWN_SUCCESS, MPIMethods::getWorldSize() - 1);
  } else {
    workerActive = false;
  }
}

void LearnerSGDEOnOffParallel::workBatch(Dataset dataset, size_t batchOffset,
                                         bool doCrossValidation) {
  waitForAllGridsConsistent();

  // assemble next batch
  std::cout << "Learning with batch of size " << dataset.getNumberInstances() << " at offset "
            << batchOffset << std::endl;
  assembleNextBatchData(&dataset, &batchOffset);
  D(std::cout << "Batch of size " << dataset.getNumberInstances()
              << " assembled, starting with training." << std::endl;)

  // train the model with current batch
  train(dataset, doCrossValidation);

  // Batch offset was already modified by assembleNextBatch
  D(std::cout << "Batch " << batchOffset - dataset.getNumberInstances() << " completed."
              << std::endl;)
  for (size_t classIndex = 0; classIndex < getNumClasses(); classIndex++) {
    D(std::cout << "Sending alpha values to master for class " << classIndex
                << " with grid version " << getLocalGridVersion(classIndex) << std::endl;)
    DataVector alphaVector = *(alphas[classIndex]);
    MPIMethods::sendMergeGridNetworkMessage(classIndex, batchOffset, dataset.getNumberInstances(),
                                            alphaVector);

    D(DataVector &dataVector = getDensityFunctions()[classIndex].first->getAlpha();
      std::cout << "Local alpha sum " << classIndex << " is now "
                << std::accumulate(dataVector.begin(), dataVector.end(), 0.0) << std::endl;)
  }
  D(std::cout << "Completed work batch " << batchOffset - dataset.getNumberInstances()
              << " requested by master." << std::endl;)
}

void LearnerSGDEOnOffParallel::waitForAllGridsConsistent() {
  size_t classIndex = 0;
  while (classIndex < localGridVersions.size()) {
    // We need to wait if the grid is not consistent or when there are differing
    // grid versions
    if (!checkGridStateConsistent(classIndex) ||
        getLocalGridVersion(classIndex) != getLocalGridVersion(0)) {
      std::cout << "Attempted to train from an inconsistent grid " << classIndex << " version "
                << getLocalGridVersion(classIndex) << std::endl;
      MPIMethods::waitForGridConsistent(classIndex);
      // start over, waiting might have changed other grids
      classIndex = 0;
    } else {
      classIndex++;
    }
  }
}

bool LearnerSGDEOnOffParallel::isVersionConsistent(size_t version) {
  return version >= MINIMUM_CONSISTENT_GRID_VERSION;
}

size_t LearnerSGDEOnOffParallel::assignBatchToWorker(size_t batchOffset, bool doCrossValidation) {
  AssignTaskResult assignTaskResult{};
  mpiTaskScheduler.assignTaskVariableTaskSize(TRAIN_FROM_BATCH, assignTaskResult);

  if (assignTaskResult.taskSize + batchOffset > trainData.getNumberInstances()) {
    std::cout << "Shortening last batch." << std::endl;
    assignTaskResult.taskSize = trainData.getNumberInstances() - batchOffset;
  }

  std::cout << "Assigning batch " << batchOffset << " to worker " << assignTaskResult.workerID
            << " with size " << assignTaskResult.taskSize << std::endl;
  MPIMethods::assignBatch(assignTaskResult.workerID, batchOffset, assignTaskResult.taskSize,
                          doCrossValidation);
  return assignTaskResult.taskSize;
}

void LearnerSGDEOnOffParallel::mergeAlphaValues(size_t classIndex, size_t remoteGridVersion,
                                                DataVector dataVector, size_t batchOffset,
                                                size_t batchSize, bool isLastPacketInSeries) {
  MPIMethods::waitForGridConsistent(classIndex);

  D(std::cout << "Remote alpha sum " << classIndex << " is "
              << std::accumulate(dataVector.begin(), dataVector.end(), 0.0) << std::endl;
    std::cout << "Batch size is " << batchSize << std::endl;)

  if (!isVersionConsistent(remoteGridVersion)) {
    std::cout << "Received merge request with inconsistent grid " << classIndex << " version "
              << remoteGridVersion << std::endl;
    throw algorithm_exception("Received a merge request to an inconsistent grid");
  }

  size_t localGridVersion = getLocalGridVersion(classIndex);
  if (isLastPacketInSeries) {
    mpiTaskScheduler.onMergeRequestIncoming(batchOffset, batchSize, remoteGridVersion,
                                            localGridVersion);
  }

  if (remoteGridVersion != localGridVersion) {
    D(std::cout << "Received merge grid request with incorrect grid version!"
                << " local: " << localGridVersion << ", remote: " << remoteGridVersion
                << std::endl;)
    if (remoteGridVersion + 1 == localGridVersion) {
      RefinementResult &refinementResult = refinementHandler.getRefinementResult(classIndex);
      std::list<size_t> &deletedPoints = refinementResult.deletedGridPointsIndices;
      std::list<LevelIndexVector> &addedPoints = refinementResult.addedGridPoints;

      D(std::cout << "Attempting to automatically compensate for outdated grid." << std::endl
                  << "Refinement result has " << addedPoints.size() << " additions and "
                  << deletedPoints.size() << " deletions" << std::endl
                  << "The original remote size is " << dataVector.size() << std::endl;)
      if (!deletedPoints.empty() || !addedPoints.empty()) {
        D(std::cout << "Found necessary refinement data" << std::endl;)

        // See DBMatOnlineDe::updateAlpha()
        if (!deletedPoints.empty()) {
          D(std::cout << "Copying vector (deleting deleted grid points)." << std::endl;)
          DataVector newAlpha(dataVector.getSize() - deletedPoints.size() + addedPoints.size());
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
        throw sgpp::base::algorithm_exception("Missing refinement data for alpha update.");
      }
    } else {
      std::cout << "Merge request " << batchOffset << ", size " << batchSize
                << ", older than one refinement cycle. Increase the refinement "
                   "period."
                << std::endl;
      throw sgpp::base::algorithm_exception("Merge request older than one refinement cycle.");
    }
  }

  DataVector &localAlpha = *(alphas[classIndex]);
  if (localAlpha.size() != dataVector.size()) {
    std::cout << "Received merge request with incorrect size (local " << localAlpha.size()
              << ", remote " << dataVector.size() << "), local version is "
              << localGridVersions[classIndex] << std::endl;
    throw sgpp::base::algorithm_exception("Merge request with incorrect size received.");
  }

  if (usePrior) {
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

size_t LearnerSGDEOnOffParallel::getLocalGridVersion(size_t classIndex) {
  return localGridVersions[classIndex];
}

bool LearnerSGDEOnOffParallel::checkAllGridsConsistent() {
  return std::all_of(localGridVersions.begin(), localGridVersions.end(),
                     [](size_t version) { return isVersionConsistent(version); });
}

bool LearnerSGDEOnOffParallel::checkGridStateConsistent(size_t classIndex) {
  if (localGridVersions.size() < classIndex) {
    throw algorithm_exception("Received request for consistency of class " + classIndex);
  }
  return isVersionConsistent(localGridVersions[classIndex]);
}

void LearnerSGDEOnOffParallel::setLocalGridVersion(size_t classIndex, size_t gridVersion) {
  D(if (localGridVersions[classIndex] != gridVersion) {
    std::cout << "Grid " << classIndex << " now has version " << gridVersion << " (previously "
              << localGridVersions[classIndex] << ")" << std::endl;
  })
  if (!checkGridStateConsistent(classIndex) && isVersionConsistent(gridVersion)) {
    std::cout << "Grid " << classIndex << " has been fully updated to version " << gridVersion
              << std::endl;
  }
  localGridVersions[classIndex] = gridVersion;
}

MPITaskScheduler &LearnerSGDEOnOffParallel::getScheduler() { return mpiTaskScheduler; }

std::unique_ptr<DBMatOffline> &LearnerSGDEOnOffParallel::getOffline() { return offline; }

Dataset &LearnerSGDEOnOffParallel::getTrainData() { return trainData; }

Dataset *LearnerSGDEOnOffParallel::getValidationData() { return validationData; }

RefinementHandler &LearnerSGDEOnOffParallel::getRefinementHandler() { return refinementHandler; }
}  // namespace datadriven
}  // namespace sgpp

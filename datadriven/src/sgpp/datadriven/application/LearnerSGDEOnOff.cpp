// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEOrthoAdapt.hpp>
#include <sgpp/datadriven/application/LearnerSGDEOnOff.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>
#include <sgpp/datadriven/algorithm/GridFactory.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitor.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorPeriodic.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorConvergence.hpp>

#include <chrono>
#include <cmath>
#include <limits>
#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>

using sgpp::base::Grid;
using sgpp::base::GridStorage;
using sgpp::base::GridGenerator;
using sgpp::base::OperationEval;
using sgpp::base::data_exception;
using sgpp::base::SurplusRefinementFunctor;
using sgpp::datadriven::GridFactory;

namespace sgpp {
namespace datadriven {

LearnerSGDEOnOff::LearnerSGDEOnOff(
    sgpp::base::RegularGridConfiguration& gridConfig,
    sgpp::base::AdpativityConfiguration& adaptivityConfig,
    sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
    Dataset& trainData,
    Dataset& testData, Dataset* validationData,
    DataVector& classLabels, size_t numClassesInit, bool usePrior,
    double beta, std::string matrixfile)
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
      avgErrors{0} {

  // initialize offline object
  if (matrixfile.empty()) {
    // Create an offline object that serves as template for all classes
    offline = std::unique_ptr<DBMatOffline>{
      DBMatOfflineFactory::buildOfflineObject(gridConfig,
            adaptivityConfig, regularizationConfig, densityEstimationConfig)
    };
  } else {
    offline = std::unique_ptr<DBMatOffline>{DBMatOfflineFactory::buildFromFile(matrixfile)};
  }

  // Create a grid TODO(fuchsgruber): Maybe remove the learner class entirely??
  grids.reserve(numClasses);
  densityFunctions.reserve(numClasses);
  offlineContainer.reserve(numClasses);
  alphas.reserve(numClasses);
  GridFactory gridFactory;
  for (size_t classIndex = 0; classIndex < numClasses; classIndex++) {
    // Create a grid
    std::unique_ptr<Grid> grid = std::unique_ptr<Grid> {
      gridFactory.createGrid(gridConfig, std::vector<std::vector <size_t>>())
    };
    std::unique_ptr<DBMatOffline> offlineCloned =
        std::unique_ptr<DBMatOffline>{offline->clone()};
    offlineCloned->buildMatrix(grid.get(), regularizationConfig);
    offlineCloned->decomposeMatrix(regularizationConfig, densityEstimationConfig);
    auto densEst = std::unique_ptr<DBMatOnlineDE>{
      DBMatOnlineDEFactory::buildDBMatOnlineDE(*offlineCloned, *grid,
          regularizationConfig.lambda_, beta)};
    densityFunctions.emplace_back(std::make_pair(std::move(densEst), classIndex));
    prior.emplace(classLabels[classIndex], 0.0);
    DataVector* alpha = new DataVector(offlineCloned->getGridSize());
    alphas.emplace_back(alpha);
    grids.emplace_back(std::move(grid));
    offlineContainer.emplace_back(std::move(offlineCloned));
  }
}


void LearnerSGDEOnOff::train(size_t batchSize, size_t maxDataPasses, std::string refType,
                             std::string refMonitor, size_t refPeriod, double accDeclineThreshold,
                             size_t accDeclineBufferSize, size_t minRefInterval, bool enableCv,
                             size_t nextCvStep,
                             sgpp::base::AdpativityConfiguration& adaptivityConfig,
                             sgpp::datadriven::DensityEstimationConfiguration&
                             densityEstimationConfig) {
  // counts total number of processed data points
  size_t totalInstances = 0;
  // contains list of removed grid points and number of added grid points
  // (is updated in each refinement/coarsening step)
  std::vector<std::pair<std::list<size_t>, size_t>> refineCoarse(numClasses);

  // initialize counter for dataset passes
  size_t cntDataPasses = 0;

  // initialize refinement variables
  double currentValidError = 0.0;
  double currentTrainError = 0.0;
  // create convergence monitor object
  RefinementMonitor *monitor = nullptr;
  if (refType == "periodic") {
    monitor = new RefinementMonitorPeriodic(refPeriod);
  } else {
    monitor = new RefinementMonitorConvergence(
      accDeclineThreshold, accDeclineBufferSize, minRefInterval);
  }

  // counts number of performed refinements
  size_t refCnt = 0;

  // coarsening
  // size_t coarseCnt = 0;
  // size_t maxCoarseNum = 5;
  // size_t coarsePeriod = 50;
  // size_t coarseNumPoints = 1;
  // double coarseThreshold = 1.0;

  auto& onlineObjects = getDensityFunctions();

  size_t dim = trainData.getDimension();
  // determine number of batches to process
  size_t numBatch = trainData.getNumberInstances() / batchSize;

  // print initial grid size
  for (auto& onlineObject : onlineObjects) {
    size_t classIndex = onlineObject.second;
    std::cout << "#Initial grid size of grid for class" << classLabels[onlineObject.second] << " : "
              << grids[classIndex]->getSize() << "\n";
  }

  // auxiliary variable for accuracy (error) measurement
  double acc = getAccuracy();
  avgErrors.append(1.0 - acc);

  // main loop which performs the training process
  while (cntDataPasses < maxDataPasses) {
    std::cout << "#batch-size: " << batchSize << "\n";
    std::cout << "#batches to process: " << numBatch << "\n";
    // data point counter - determines offset when selecting next batch
    size_t cnt = 0;
    // iterate over total number of batches
    for (size_t step = 1; step <= numBatch; step++) {
      std::cout << "#processing batch: " << step << "\n";

      auto begin = std::chrono::high_resolution_clock::now();

      // check if cross-validation should be performed
      bool doCv = false;
      if (enableCv) {
        if (nextCvStep == step) {
          doCv = true;
          nextCvStep *= 5;
        }
      }
      // assemble next batch
      Dataset currentBatch(batchSize, dim);
      for (size_t j = 0; j < batchSize; j++) {
        DataVector x(dim);
        trainData.getData().getRow(j + cnt, x);
        double y = trainData.getTargets().get(j + cnt);
        currentBatch.getData().setRow(j, x);
        currentBatch.getTargets().set(j, y);
      }
      // curPair = std::pair<DataMatrix*, DataVector*>(batch, batchLabels);
      cnt += batchSize;

      // train the model with current batch
      train(currentBatch, adaptivityConfig, densityEstimationConfig, doCv, &refineCoarse);

      totalInstances += currentBatch.getNumberInstances();

      // access DBMatOnlineDE-objects of all classes in order
      // to apply adaptivity to the specific sparse grids later on

      size_t refinementsNeccessary = 0;
      if (offline->isRefineable() && (refCnt < adaptivityConfig.numRefinements_) && monitor) {
        monitor->pushToBuffer(currentBatch.getNumberInstances(),
            currentValidError, currentTrainError);
        refinementsNeccessary = monitor->refinementsNeccessary();
      }

      // if the Cholesky decomposition is chosen as factorization method
      // refinement
      // and coarsening methods can be applied
      while (refinementsNeccessary > 0) {
        std::cout << "refinement at iteration: " << totalInstances << "\n";
        refine(*monitor, adaptivityConfig, densityEstimationConfig, refineCoarse, refType);
        refCnt += 1;
        refinementsNeccessary--;
      }

      // save current error
      if (totalInstances % 10 == 0) {
        acc = getAccuracy();
        avgErrors.append(1.0 - acc);
      }

      auto end = std::chrono::high_resolution_clock::now();
      std::cout << "Processing batch in "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
                << "ms" << std::endl;
    }
    cntDataPasses++;
    processedPoints = 0;
  }  // end while

  std::cout << "#Training finished"
            << "\n";
}

void LearnerSGDEOnOff::train(
    Dataset& dataset,
    sgpp::base::AdpativityConfiguration& adaptivityConfig,
    sgpp::datadriven::DensityEstimationConfiguration& densityEstimationConfig,
    bool doCv,
    std::vector<std::pair<std::list<size_t>, size_t>>* refineCoarse) {
  size_t dim = dataset.getDimension();

  // create an empty matrix for every class:
  std::vector<std::unique_ptr<DataMatrix>> classData;
  classData.reserve(classLabels.getSize());
  std::vector<std::pair<DataMatrix*, double>> trainDataClasses;
  trainDataClasses.reserve(classLabels.getSize());

  std::map<double, int> classIndices;  // maps class labels to indices
  int index = 0;
  for (size_t i = 0; i < classLabels.getSize(); i++) {
    std::unique_ptr<DataMatrix> dataMatrix = std::make_unique<DataMatrix>(0, dim);
    auto trainDataClassPair = std::make_pair(&(*dataMatrix), classLabels[i]);
    classData.emplace_back(std::move(dataMatrix));
    trainDataClasses.emplace_back(trainDataClassPair);
    classIndices.emplace(std::make_pair(classLabels[i], index));
    index++;
  }
  // split the data into the different classes:
  for (size_t i = 0; i < dataset.getNumberInstances(); i++) {
    double classLabel = dataset.getTargets()[i];
    DataVector vec(dim);
    dataset.getData().getRow(i, vec);
    auto& p = trainDataClasses[classIndices[classLabel]];
    p.first->appendRow(vec);
  }

  // compute density functions
  train(trainDataClasses, densityEstimationConfig, doCv, refineCoarse);
}

void LearnerSGDEOnOff::train(std::vector<std::pair<DataMatrix*, double>>& trainDataClasses,
                             sgpp::datadriven::DensityEstimationConfiguration&
                            densityEstimationConfig,
                             bool doCv,
                             std::vector<std::pair<std::list<size_t>, size_t>>* refineCoarse) {
  size_t numberOfDataPoints = 0;
  for (size_t i = 0; i < trainDataClasses.size(); i++) {
    numberOfDataPoints += trainDataClasses[i].first->getSize();
  }
  for (size_t i = 0; i < trainDataClasses.size(); i++) {
    auto& p = trainDataClasses[i];

    if ((*p.first).getNrows() > 0) {
      // update density function for current class
      size_t classIndex = densityFunctions[i].second;
      densityFunctions[i].first->computeDensityFunction(
          *(alphas[i]), *p.first, *(grids[classIndex]), densityEstimationConfig, true, doCv,
          &(*refineCoarse)[i].first, (*refineCoarse)[i].second);
      (*refineCoarse)[i].first.clear();
      (*refineCoarse)[i].second = 0;

      if (usePrior) {
        this->prior[p.second] =
            ((this->prior[p.second] * static_cast<double>(processedPoints)) +
             static_cast<double>(p.first->getSize())) /
            (static_cast<double>(numberOfDataPoints) + static_cast<double>(processedPoints));
      } else {
        this->prior[p.second] = 1.;
      }
    }
  }

  this->processedPoints += numberOfDataPoints;
  trained = true;
}


double LearnerSGDEOnOff::getAccuracy() const {
  DataVector computedLabels{testData.getNumberInstances()};
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

void LearnerSGDEOnOff::predict(DataMatrix& data, DataVector& result) const {
  // calculate per class densities
  std::vector<DataVector> perClassDensities;
  for (auto& densityFunction : densityFunctions) {
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

int LearnerSGDEOnOff::predict(DataVector& p) const {
  DataMatrix ptmp(1, p.getSize());
  ptmp.setRow(0, p);
  DataVector r{ptmp.getNrows()};
  this->predict(ptmp, r);
  return static_cast<int>(r[0]);
}

double LearnerSGDEOnOff::getError(Dataset& dataset) const {
  double res = -1.0;

  DataVector computedLabels{dataset.getNumberInstances()};
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

void LearnerSGDEOnOff::storeResults() {
  DataVector classesComputed{testData.getNumberInstances()};
  predict(testData.getData(), classesComputed);

  std::ofstream output;
  // write predicted classes to csv file
  output.open("SGDEOnOff_predicted_classes.csv");
  if (output.fail()) {
    std::cout << "failed to create csv file!"
              << "\n";
  } else {
    for (size_t i = 0; i < classesComputed.getSize(); i++) {
      DataVector x(2);
      testData.getData().getRow(i, x);
      output << x[0] << ";" << x[1] << ";" << classesComputed[i] << "\n";
    }
    output.close();
  }
  // write grids to csv file
  for (size_t i = 0; i < numClasses; i++) {
    // auto densEst = densityFunctions[i].first.get(); // seems to be unused?
    size_t classIndex = densityFunctions[i].second;
    output.open("SGDEOnOff_grid_" + std::to_string(densityFunctions[i].second) + ".csv");
    if (output.fail()) {
      std::cout << "failed to create csv file!"
                << "\n";
    } else {
      GridStorage& storage = grids[classIndex]->getStorage();
      for (auto& iter : storage) {
        DataVector gpCoord(trainData.getDimension());
        storage.getCoordinates(*(iter.first), gpCoord);
        for (size_t d = 0; d < gpCoord.getSize(); d++) {
          if (d < gpCoord.getSize() - 1) {
            output << gpCoord[d] << ";";
          } else {
            output << gpCoord[d] << "\n";
          }
        }
      }
      output.close();
    }
  }
  // write density function evaluations to csv file
  double stepSize = 0.01;
  DataMatrix values(0, 2);
  DataVector range(101);
  for (size_t i = 0; i < 101; i++) {
    range.set(i, stepSize * (static_cast<double>(i)));
  }
  for (size_t i = 0; i < range.getSize(); i++) {
    for (size_t j = 0; j < range.getSize(); j++) {
      DataVector row(2);
      row.set(1, range.get(i));
      row.set(0, range.get(j));
      values.appendRow(row);
    }
  }
  // evaluate each density function at all points from values
  // and write result to csv file
  for (size_t j = 0; j < densityFunctions.size(); j++) {
    auto& pair = densityFunctions[j];
    size_t classIndex = pair.second;
    output.open("SGDEOnOff_density_fun_" + std::to_string(pair.second) + "_evals.csv");
    for (size_t i = 0; i < values.getNrows(); i++) {
      // get next test sample x
      DataVector x(2);
      values.getRow(i, x);
      double density = pair.first->eval(*(alphas[classIndex]), x, *(grids[classIndex]), true) *
          this->prior[classLabels[classIndex]];
      output << density << ";"
             << "\n";
    }
    output.close();
  }
}

void LearnerSGDEOnOff::getDensities(DataVector& point, DataVector& density) const {
  for (size_t i = 0; i < densityFunctions.size(); i++) {
    auto& pair = densityFunctions[i];
    size_t classIndex = densityFunctions[i].second;
    density[i] = pair.first->eval(*(alphas[classIndex]), point, *(grids[classIndex]));
  }
}

void LearnerSGDEOnOff::setCrossValidationParameters(int lambdaStep, double lambdaStart,
                                                    double lambdaEnd, DataMatrix* test,
                                                    DataMatrix* testRes, bool logscale) {
  for (auto& destFunction : densityFunctions) {
    destFunction.first->setCrossValidationParameters(lambdaStep, lambdaStart, lambdaEnd, test,
                                                     testRes, logscale);
  }
}

/*double LearnerSGDEOnOff::getBestLambda() {
  // return online->getBestLambda();
  return 0.; // TODO
}*/

size_t LearnerSGDEOnOff::getNumClasses() const { return numClasses; }

void LearnerSGDEOnOff::getAvgErrors(DataVector& result) const { result = avgErrors; }

ClassDensityConntainer& LearnerSGDEOnOff::getDensityFunctions() { return densityFunctions; }

void LearnerSGDEOnOff::updateAlpha(size_t classIndex, std::list<size_t>* deletedPoints,
    size_t newPoints) {
  DataVector& alpha = *(alphas[classIndex]);
  if (alpha.getSize() != 0 && deletedPoints != nullptr && !deletedPoints->empty()) {
    std::vector<size_t> deletedPoints_{std::begin(*deletedPoints), std::end(*deletedPoints)};
    DataVector newAlpha{alpha.getSize() - deletedPoints->size() + newPoints};
    for (size_t i = 0; i < alpha.getSize(); i++) {
      if (std::find(deletedPoints_.begin(), deletedPoints_.end(), i) != deletedPoints_.end()) {
        continue;
      } else {
        newAlpha.append(alpha.get(i));
      }
    }
    // set new alpha
    alpha = std::move(newAlpha);
  }
  if (newPoints > 0) {
    alpha.resizeZero(alpha.getSize() + newPoints);
  }
}

void LearnerSGDEOnOff::refine(RefinementMonitor& monitor,
                              sgpp::base::AdpativityConfiguration& adaptivityConfig,
                              sgpp::datadriven::DensityEstimationConfiguration&
                              densityEstimationConfig,
                              std::vector<std::pair<std::list<size_t>, size_t>>& refineCoarse,
                              std::string& refType) {
  auto& onlineObjects = getDensityFunctions();
  DataVector* alphaWork;  // required for surplus refinement
  // auxiliary variables
  DataVector p(trainData.getDimension());

  size_t newPoints = 0;
  std::list<size_t> deletedGridPoints;

  // acc = getAccuracy();
  // avgErrors.append(1.0 - acc);

  // bundle grids and surplus vector pointer needed for refinement
  // (for zero-crossings refinement, data-based refinement)
  std::vector<Grid*> gridsVector;
  bool levelPenalize = false;  // Multiplies penalzing term for fine levels
  bool preCompute = true;      // Precomputes and caches evals for zrcr
  MultiGridRefinementFunctor* func = nullptr;
  gridsVector.reserve(grids.size());
  for (size_t idx = 0; idx < grids.size(); idx++) {
    gridsVector.push_back(&(*(grids[idx])));
  }

  // Zero-crossing-based refinement
  ZeroCrossingRefinementFunctor funcZrcr{gridsVector, alphas, adaptivityConfig.noPoints_,
                                         levelPenalize, preCompute};

  // Data-based refinement. Needs a problem dependent coeffA. The values
  // can be determined by testing (aim at ~10 % of the training data is
  // to be marked relevant). Cross-validation or similar can/should be
  // employed
  // to determine this value.
  std::vector<double> coeffA;
  coeffA.push_back(1.2);  // ripley 1.2
  coeffA.push_back(1.2);  // ripley 1.2
  DataMatrix* trainDataRef = &(trainData.getData());
  DataVector* trainLabelsRef = &(trainData.getTargets());
  DataBasedRefinementFunctor funcData = DataBasedRefinementFunctor{
      gridsVector, alphas, trainDataRef, trainLabelsRef, adaptivityConfig.noPoints_,
      levelPenalize, coeffA};

  if (refType == "zero") {
    func = &funcZrcr;
  } else if (refType == "data") {
    func = &funcData;
  }

  // note: if you dont want to coarsen, just set coarseCnt and maxCoarseNum both to zero
  size_t coarseCnt = 0;
  size_t maxCoarseNum = 3;
  // size_t coarsePeriod = 50;
  size_t coarseNumPoints = 15;
  double coarseThreshold = 1.0;

  if (offline->isRefineable()) {
    // perform refinement/coarsening for each grid
    for (size_t idx = 0; idx < getNumClasses(); idx++) {
      // perform refinement/coarsening for grid which corresponds to current
      // index
      std::cout << "\nRefinement and coarsening for class: " << idx << "\n";
      auto densEst = onlineObjects[idx].first.get();
      std::cout << "Size before adaptivity: " << grids[idx]->getSize() << "\n";

      GridGenerator& gridGen = grids[idx]->getGenerator();

      size_t sizeBeforeRefine = grids[idx]->getSize();
      size_t sizeAfterRefine = grids[idx]->getSize();

      if (refType == "surplus") {
        std::unique_ptr<OperationEval> opEval(op_factory::createOperationEval(*(grids[idx])));
        GridStorage& gridStorage = grids[idx]->getStorage();
        alphaWork = alphas[idx];
        DataVector alphaWeight(alphaWork->getSize());
        // determine surpluses
        for (size_t k = 0; k < gridStorage.getSize(); k++) {
          // sets values of p to the coordinates of the given GridPoint gp
          gridStorage.getPoint(k).getStandardCoordinates(p);
          // multiply k-th alpha with the evaluated function at grind-point
          // k
          alphaWeight[k] = alphaWork->get(k) * opEval->eval(*alphaWork, p);
        }

        // ### begin: comment this out if a non-surplus based strategy is used
        // Perform Coarsening (surplus based)

        // forbid coarsening of initial gridpoints in case of OrthoAdapt
        size_t minIndexAllowed = 0;
        if (densityEstimationConfig.decomposition_ ==
            sgpp::datadriven::MatrixDecompositionType::OrthoAdapt) {
          minIndexAllowed =
              static_cast<sgpp::datadriven::DBMatOfflineOrthoAdapt&>(*offline).getGridSize();
        }

        if (coarseCnt < maxCoarseNum) {
          sgpp::base::HashCoarsening coarse_;
          // std::cout << "\n" << "Start coarsening\n";

          // Coarsening based on surpluses
          sgpp::base::SurplusCoarseningFunctor scf(alphaWeight, coarseNumPoints, coarseThreshold);

          std::cout << "Size before coarsening:" << grids[idx]->getSize() << "\n";

          // int old_size = grid->getSize();

          // std::cout << "minIndexAllowed: " << minIndexAllowed << std::endl;
          std::cout << "gridSize: " << grids[idx]->getSize() << std::endl;
          coarse_.free_coarsen_NFirstOnly(grids[idx]->getStorage(), scf, alphaWeight,
              grids[idx]->getSize(), minIndexAllowed);

          std::cout << "Size after coarsening:" << grids[idx]->getSize() << "\n";
          // int new_size = grid->getSize();

          deletedGridPoints.clear();
          deletedGridPoints = coarse_.getDeletedPoints();

          (refineCoarse)[idx].first = deletedGridPoints;

          coarseCnt++;
        }

        // perform refinement (surplus based)
        sizeBeforeRefine = grids[idx]->getSize();
        std::cout << "Size before refine: " << sizeBeforeRefine << std::endl;
        // simple refinement based on surpluses
        SurplusRefinementFunctor srf(alphaWeight, adaptivityConfig.noPoints_);
        if (offline->interactions.size() == 0) {
          gridGen.refine(srf);
        } else {
          gridGen.refineInter(srf, offline->interactions);
        }
        sizeAfterRefine = grids[idx]->getSize();
      } else if ((refType == "data") || (refType == "zero")) {
        if (preCompute) {
          // precompute the evals (needs to be done once per step, before
          // any refinement is done
          func->preComputeEvaluations();
        }
        func->setGridIndex(idx);
        // perform refinement (zero-crossings-based / data-based)
        sizeBeforeRefine = grids[idx]->getSize();
        if (offline->interactions.size() == 0) {
          gridGen.refine(*func);
        } else {
          gridGen.refineInter(*func, offline->interactions);
        }
        sizeAfterRefine = grids[idx]->getSize();
      }

      std::cout << "grid size after refine: " << grids[idx]->getSize() << "\n";

      newPoints = sizeAfterRefine - sizeBeforeRefine;
      // std::cout << "will be adding " << newPoints << " new points\n";
      refineCoarse[idx].second = newPoints;

      // TODO(Kilian) the .updateSystemMatrixDecomposition function
      // of online_ortho_adapt is currently designed to
      // return a list of gridpoints which weren't allowed to be coarsened. But it seems
      // appropriate to redesign the functors in a way, that already considers these points
      // when coarsening the grid itself.
      densEst->updateSystemMatrixDecomposition(densityEstimationConfig,
          *(grids[idx]), newPoints, deletedGridPoints, densEst->getBestLambda());


      // update alpha vector
      updateAlpha(idx, &refineCoarse[idx].first, refineCoarse[idx].second);
    }
  }
}

}  // namespace datadriven
}  // namespace sgpp


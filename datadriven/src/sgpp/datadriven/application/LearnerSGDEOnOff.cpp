// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

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
#include <sgpp/datadriven/application/LearnerSGDEOnOff.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>

#include <cmath>
#include <ctime>
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

namespace sgpp {
namespace datadriven {

LearnerSGDEOnOff::LearnerSGDEOnOff(DBMatDensityConfiguration& dconf, Dataset& trainData,
                                   Dataset& testData, Dataset* validationData,
                                   DataVector& classLabels, size_t numClassesInit, bool usePrior,
                                   double beta, double lambda)
    : trainData{trainData},
      testData{testData},
      validationData{validationData},
      classLabels{classLabels},
      numClasses{numClassesInit},
      usePrior{usePrior},
      prior{},
      beta{beta},
      trained{false},
      offline{nullptr},
      offlineContainer{},
      densityFunctions{},
      processedPoints{0},
      avgErrors{0} {
  // initialize offline object
  offline = std::unique_ptr<DBMatOffline>{DBMatOfflineFactory::buildOfflineObject(dconf)};
  offline->buildMatrix();
  offline->decomposeMatrix();

  // initialize density functions for each class
  densityFunctions.reserve(numClasses);
  // if the Cholesky decomposition is chosen declare separate Online-objects for
  // every class
  if (offline->getConfig().decomp_type_ == DBMatDecompostionType::DBMatDecompChol) {
    offlineContainer.reserve(numClasses);
    // every class gets his own online object
    for (size_t classIndex = 0; classIndex < numClasses; classIndex++) {
      offlineContainer.emplace_back(
          std::make_unique<DBMatOfflineChol>(static_cast<DBMatOfflineChol&>(*offline)));
      auto densEst = std::make_unique<DBMatOnlineDE>(*(offlineContainer.back()), beta);
      densityFunctions.emplace_back(std::make_pair(std::move(densEst), classLabels[classIndex]));
    }
  } else {
    densityFunctions.reserve(numClasses);
    for (size_t classIndex = 0; classIndex < numClasses; classIndex++) {
      auto densEst = std::make_unique<DBMatOnlineDE>(*offline, beta);
      densityFunctions.emplace_back(std::make_pair(std::move(densEst), classLabels[classIndex]));
    }
  }

  for (size_t i = 0; i < numClasses; i++) {
    prior.emplace(classLabels[i], 0.0);
  }

  for (auto& iter : densityFunctions) {
    iter.first->setLambda(lambda);
  }
}

void LearnerSGDEOnOff::train(size_t batchSize, size_t maxDataPasses, std::string refType,
                             std::string refMonitor, size_t refPeriod, double accDeclineThreshold,
                             size_t accDeclineBufferSize, size_t minRefInterval, bool enableCv,
                             size_t nextCvStep) {
  // counts total number of processed data points
  size_t totalInstances = 0;
  // contains list of removed grid points and number of added grid points
  // (is updated in each refinement/coarsening step)
  std::vector<std::pair<std::list<size_t>, size_t> >* refineCoarse =
      new std::vector<std::pair<std::list<size_t>, size_t> >(numClasses);

  // auxiliary variables
  DataVector* alphaWork;  // required for surplus refinement
  DataVector p(trainData.getDimension());

  // initialize counter for dataset passes
  size_t cntDataPasses = 0;

  // initialize refinement variables
  double currentValidError = 0.0;
  double currentTrainError = 0.0;
  // create convergence monitor object
  std::shared_ptr<ConvergenceMonitor> monitor(
      new ConvergenceMonitor(accDeclineThreshold, accDeclineBufferSize, minRefInterval));
  bool doRefine = false;  // is set to 'true' by refinement monitor
  // counts number of performed refinements
  size_t refCnt = 0;

  // coarsening
  // size_t coarseCnt = 0;
  // size_t maxCoarseNum = 5;
  // size_t coarsePeriod = 50;
  // size_t coarseNumPoints = 1;
  // double coarseThreshold = 1.0;

  std::list<size_t> deletedGridPoints;
  size_t newPoints = 0;

  auto& onlineObjects = getDensityFunctions();

  size_t dim = trainData.getDimension();
  // determine number of batches to process
  size_t numBatch = trainData.getNumberInstances() / batchSize;

  // print initial grid size
  for (auto& onlineObject : onlineObjects) {
    auto densEst = onlineObject.first.get();
    Grid& grid = densEst->getOfflineObject().getGrid();
    std::cout << "#Initial grid size of grid for class" << onlineObject.second << " : "
              << grid.getSize() << "\n";
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
      train(currentBatch, doCv, refineCoarse);

      totalInstances += currentBatch.getNumberInstances();

      // access DBMatOnlineDE-objects of all classes in order
      // to apply adaptivity to the specific sparse grids later on

      // check if refinement should be performed
      if (refMonitor == "periodic") {
        // check periodic monitor
        if ((offline->getConfig().decomp_type_ == DBMatDecompostionType::DBMatDecompChol) &&
            (totalInstances > 0) && (totalInstances % refPeriod == 0) &&
            (refCnt < offline->getConfig().numRefinements_)) {
          doRefine = true;
        }
      } else if (refMonitor == "convergence") {
        // check convergence monitor
        if (validationData == nullptr) {
          throw base::data_exception("No validation data for checking convergence provided!");
        }
        if ((offline->getConfig().decomp_type_ == DBMatDecompostionType::DBMatDecompChol) &&
            (refCnt < offline->getConfig().numRefinements_)) {
          currentValidError = getError(*validationData);
          currentTrainError = getError(trainData);  // if train dataset is large
                                                    // use a subset for error
                                                    // evaluation
          monitor->pushToBuffer(currentValidError, currentTrainError);
          if (monitor->nextRefCnt > 0) {
            monitor->nextRefCnt--;
          }
          if (monitor->nextRefCnt == 0) {
            doRefine = monitor->checkConvergence();
          }
        }
      }
      // if the Cholesky decomposition is chosen as factorization method
      // refinement
      // and coarsening methods can be applied
      if (doRefine) {
        // acc = getAccuracy();
        // avgErrors.append(1.0 - acc);
        std::cout << "refinement at iteration: " << totalInstances << "\n";
        // bundle grids and surplus vector pointer needed for refinement
        // (for zero-crossings refinement, data-based refinement)
        std::vector<Grid*> grids;
        std::vector<DataVector*> alphas;
        for (size_t i = 0; i < getNumClasses(); i++) {
          auto densEst = onlineObjects[i].first.get();
          grids.push_back(&(densEst->getOfflineObject().getGrid()));
          alphas.push_back(densEst->getAlpha());
        }
        bool levelPenalize = false;  // Multiplies penalzing term for fine levels
        bool preCompute = true;      // Precomputes and caches evals for zrcr
        sgpp::datadriven::MultiGridRefinementFunctor* func = nullptr;

        // Zero-crossing-based refinement
        sgpp::datadriven::ZeroCrossingRefinementFunctor funcZrcr =
            *(new sgpp::datadriven::ZeroCrossingRefinementFunctor(
                grids, alphas, offline->getConfig().ref_noPoints_, levelPenalize, preCompute));

        // Data-based refinement. Needs a problem dependent coeffA. The values
        // can be determined by testing (aim at ~10 % of the training data is
        // to be marked relevant). Cross-validation or similar can/should be
        // employed
        // to determine this value.
        std::vector<double> coeffA;
        coeffA.push_back(1.2);  // ripley 1.2
        coeffA.push_back(1.2);  // ripley 1.2
        base::DataMatrix* trainDataRef = &(trainData.getData());
        base::DataVector* trainLabelsRef = &(trainData.getTargets());
        sgpp::datadriven::DataBasedRefinementFunctor funcData =
            *(new sgpp::datadriven::DataBasedRefinementFunctor(
                grids, alphas, trainDataRef, trainLabelsRef, offline->getConfig().ref_noPoints_,
                levelPenalize, coeffA));
        if (refType == "zero") {
          func = &funcZrcr;
        } else if (refType == "data") {
          func = &funcData;
        }

        // perform refinement/coarsening for each grid
        for (size_t idx = 0; idx < getNumClasses(); idx++) {
          // perform refinement/coarsening for grid which corresponds to current
          // index
          std::cout << "Refinement and coarsening for class: " << idx << "\n";
          auto densEst = onlineObjects[idx].first.get();
          Grid& grid = densEst->getOfflineObject().getGrid();
          std::cout << "Size before adaptivity: " << grid.getSize() << "\n";

          GridGenerator& gridGen = grid.getGenerator();

          size_t sizeBeforeRefine = grid.getSize();
          size_t sizeAfterRefine = grid.getSize();

          if (refType == "surplus") {
            std::unique_ptr<OperationEval> opEval(op_factory::createOperationEval(grid));
            GridStorage& gridStorage = grid.getStorage();
            alphaWork = densEst->getAlpha();
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
            sizeBeforeRefine = grid.getSize();
            // simple refinement based on surpluses
            SurplusRefinementFunctor srf(alphaWeight, offline->getConfig().ref_noPoints_);
            gridGen.refine(srf);
            sizeAfterRefine = grid.getSize();
          } else if ((refType == "data") || (refType == "zero")) {
            if (preCompute) {
              // precompute the evals (needs to be done once per step, before
              // any refinement is done
              func->preComputeEvaluations();
            }
            func->setGridIndex(idx);
            // perform refinement (zero-crossings-based / data-based)
            sizeBeforeRefine = grid.getSize();
            gridGen.refine(*func);
            sizeAfterRefine = grid.getSize();
          }

          std::cout << "grid size after adaptivity: " << grid.getSize() << "\n";

          newPoints = sizeAfterRefine - sizeBeforeRefine;
          (*refineCoarse)[idx].second = newPoints;
          // apply grid changes to the Cholesky factorization
          if (offline->getConfig().decomp_type_ == DBMatDecompostionType::DBMatDecompChol) {
            static_cast<DBMatOfflineChol&>(densEst->getOfflineObject())
                .choleskyModification(newPoints, deletedGridPoints, densEst->getBestLambda());
          }
          // update alpha vector
          densEst->updateAlpha(&(*refineCoarse)[idx].first, (*refineCoarse)[idx].second);
        }
        refCnt += 1;
        doRefine = false;
        if (refMonitor == "convergence") {
          monitor->nextRefCnt = monitor->minRefInterval;
        }
      }

      // save current error
      if (totalInstances % 10 == 0) {
        acc = getAccuracy();
        avgErrors.append(1.0 - acc);
      }
    }
    cntDataPasses++;
    processedPoints = 0;
  }  // end while

  std::cout << "#Training finished"
            << "\n";

  // delete offline;
  delete refineCoarse;
}

void LearnerSGDEOnOff::train(Dataset& dataset, bool doCv,
                             std::vector<std::pair<std::list<size_t>, size_t> >* refineCoarse) {
  size_t dim = dataset.getDimension();

  // create an empty matrix for every class:
  std::vector<std::pair<DataMatrix*, double> > trainDataClasses;
  trainDataClasses.reserve(classLabels.getSize());

  std::map<double, int> classIndices;  // maps class labels to indices
  int index = 0;
  for (size_t i = 0; i < classLabels.getSize(); i++) {
    DataMatrix* m = new DataMatrix(0, dim);
    trainDataClasses.emplace_back(std::make_pair(m, classLabels[i]));
    classIndices.emplace(std::make_pair(classLabels[i], index));
    index++;
  }
  // split the data into the different classes:
  for (size_t i = 0; i < dataset.getNumberInstances(); i++) {
    double classLabel = dataset.getTargets()[i];
    DataVector vec(dim);
    dataset.getData().getRow(i, vec);
    std::pair<DataMatrix*, double> p = trainDataClasses[classIndices[classLabel]];
    p.first->appendRow(vec);
  }
  // compute density functions
  train(trainDataClasses, doCv, refineCoarse);

  // delete DataMatrix pointers:
  for (size_t i = 0; i < trainDataClasses.size(); i++) {
    delete trainDataClasses[i].first;
  }
}

void LearnerSGDEOnOff::train(std::vector<std::pair<DataMatrix*, double> >& trainDataClasses,
                             bool doCv,
                             std::vector<std::pair<std::list<size_t>, size_t> >* refineCoarse) {
  size_t numberOfDataPoints = 0;
  for (size_t i = 0; i < trainDataClasses.size(); i++) {
    numberOfDataPoints += trainDataClasses[i].first->getSize();
  }
  for (size_t i = 0; i < trainDataClasses.size(); i++) {
    std::pair<DataMatrix*, double> p = trainDataClasses[i];

    if ((*p.first).getNrows() > 0) {
      // update density function for current class
      densityFunctions[i].first->computeDensityFunction(
          *p.first, true, doCv, &(*refineCoarse)[i].first, (*refineCoarse)[i].second);
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
  DataVector computedLabels = predict(testData.getData());
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

base::DataVector LearnerSGDEOnOff::predict(DataMatrix& data) const {
  base::DataVector result(data.getNrows());

  /*if(not trained) {
    std::cerr << "LearnerSGDEOnOff: Not trained!\n";
    exit(-1);
  }*/

  for (size_t i = 0; i < data.getNrows(); i++) {
    double max = std::numeric_limits<double>::max() * (-1);
    double max_class = 0;
    // compute the maximum density:
    DataVector p(data.getNcols());
    data.getRow(i, p);
    for (size_t j = 0; j < densityFunctions.size(); j++) {
      auto& pair = densityFunctions[j];
      // double density = pair.first->eval(p)*this->prior[pair.second];
      double density = pair.first->eval(p, true) * prior.at(pair.second);
      if (density > max) {
        max = density;
        max_class = pair.second;
      }
    }
    if (max_class == 0) {
      std::cerr << "LearnerSGDEOnOff: Warning: no best class found!"
                << "\n";
    }
    result[i] = max_class;
  }
  return result;
}

int LearnerSGDEOnOff::predict(DataVector& p) const {
  DataMatrix ptmp(1, p.getSize());
  ptmp.setRow(0, p);
  DataVector r = this->predict(ptmp);
  return static_cast<int>(r[0]);
}

double LearnerSGDEOnOff::getError(Dataset& dataset) const {
  double res = -1.0;

  DataVector computedLabels = predict(dataset.getData());
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
  DataVector classesComputed = predict(testData.getData());

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
    auto densEst = densityFunctions[i].first.get();
    Grid& grid = densEst->getOfflineObject().getGrid();
    output.open("SGDEOnOff_grid_" + std::to_string(densityFunctions[i].second) + ".csv");
    if (output.fail()) {
      std::cout << "failed to create csv file!"
                << "\n";
    } else {
      GridStorage& storage = grid.getStorage();
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
    output.open("SGDEOnOff_density_fun_" + std::to_string(pair.second) + "_evals.csv");
    for (size_t i = 0; i < values.getNrows(); i++) {
      // get next test sample x
      DataVector x(2);
      values.getRow(i, x);
      double density = pair.first->eval(x, true) * this->prior[pair.second];
      output << density << ";"
             << "\n";
    }
    output.close();
  }
}

DataVector LearnerSGDEOnOff::getDensities(DataVector& point) const {
  base::DataVector result(densityFunctions.size());
  for (size_t i = 0; i < densityFunctions.size(); i++) {
    auto& pair = densityFunctions[i];
    result[i] = pair.first->eval(point);
  }
  return result;
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

}  // namespace datadriven
}  // namespace sgpp

#endif /* USE_GSL */

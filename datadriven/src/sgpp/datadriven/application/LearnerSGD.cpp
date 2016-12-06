// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/functors/ImpurityRefinementIndicator.hpp>
#include <sgpp/base/grid/generation/functors/PredictiveRefinementIndicator.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/ImpurityRefinement.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/PredictiveRefinement.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/algorithm/ConvergenceMonitor.hpp>
#include <sgpp/datadriven/application/LearnerSGD.hpp>

#include <cmath>
#include <string>
#include <algorithm>

using sgpp::base::GridStorage;
using sgpp::base::HashRefinement;
using sgpp::base::PredictiveRefinement;
using sgpp::base::PredictiveRefinementIndicator;
using sgpp::base::ImpurityRefinement;
using sgpp::base::ImpurityRefinementIndicator;

namespace sgpp {
namespace datadriven {

LearnerSGD::LearnerSGD(sgpp::base::RegularGridConfiguration& gridConfig,
                       sgpp::base::AdpativityConfiguration& adaptivityConfig)
    : grid(nullptr),
      alpha(nullptr),
      alphaAvg(nullptr),
      trainData(nullptr),
      trainLabels(nullptr),
      testData(nullptr),
      testLabels(nullptr),
      batchData(nullptr),
      batchLabels(nullptr),
      batchError(nullptr),
      gridConfig(gridConfig),
      adaptivityConfig(adaptivityConfig),
      lambda(0),
      gamma(0),
      currentGamma(0),
      batchSize(0) {}

LearnerSGD::~LearnerSGD() {}

void LearnerSGD::initialize(sgpp::base::DataMatrix& pTrainData,
                            sgpp::base::DataVector& pTrainLabels,
                            sgpp::base::DataMatrix& pTestData,
                            sgpp::base::DataVector& pTestLabels,
                            std::shared_ptr<sgpp::base::DataMatrix> pValData,
                            std::shared_ptr<sgpp::base::DataVector> pValLabels,
                            double lambda, double gamma, size_t batchSize,
                            bool useValidData) {
  trainData = std::make_shared<base::DataMatrix>(pTrainData);
  trainLabels = std::make_shared<base::DataVector>(pTrainLabels);
  testData = std::make_shared<base::DataMatrix>(pTestData);
  testLabels = std::make_shared<base::DataVector>(pTestLabels);

  this->lambda = lambda;
  this->gamma = gamma;
  currentGamma = gamma;

  this->batchSize = batchSize;

  this->useValidData = useValidData;

  // if no validation data is provided -> create buffer
  // which contains already processed data points
  // (required for computing error contributions used for predictive refinement)
  if (!useValidData) {
    batchData = std::shared_ptr<base::DataMatrix>(
        new base::DataMatrix(0, trainData->getNcols()));
    batchData->addSize(batchSize);
    batchData->setAll(0.0);
    batchLabels =
        std::shared_ptr<base::DataVector>(new base::DataVector(batchSize));
    batchLabels->setAll(0.0);
  } else {
    batchData = pValData;
    batchLabels = pValLabels;
  }

  // vector containing error contributions for predictive refinement
  batchError =
      std::shared_ptr<base::DataVector>(new base::DataVector(batchSize));
  batchError->setAll(0.0);

  // create sparse grid
  gridConfig.dim_ = trainData->getNcols();
  grid = createRegularGrid();
  std::cout << "# initial grid size: " << grid->getSize() << std::endl;
  // surplus vector
  alpha =
      std::shared_ptr<base::DataVector>(new base::DataVector(grid->getSize()));
  alpha->setAll(0.0);
  // vector for averaged surpluses
  alphaAvg =
      std::shared_ptr<base::DataVector>(new base::DataVector(grid->getSize()));
  alphaAvg->setAll(0.0);
}

std::shared_ptr<base::Grid> LearnerSGD::createRegularGrid() {
  // load grid
  std::unique_ptr<base::Grid> uGrid;
  if (gridConfig.type_ == base::GridType::Linear) {
    uGrid.reset(base::Grid::createLinearGrid(gridConfig.dim_));
  } else if (gridConfig.type_ == base::GridType::ModLinear) {
    uGrid.reset(base::Grid::createModLinearGrid(gridConfig.dim_));
  } else {
    throw base::application_exception(
        "LearnerSGD::initialize : grid type is not supported");
  }
  uGrid->getGenerator().regular(gridConfig.level_);

  // move the grid to be shared
  std::shared_ptr<base::Grid> sGrid{std::move(uGrid)};

  return sGrid;
}

void LearnerSGD::train(size_t maxDataPasses, std::string refType,
                       std::string refMonitor, size_t refPeriod,
                       double errorDeclineThreshold,
                       size_t errorDeclineBufferSize, size_t minRefInterval) {
  size_t dim = trainData->getNcols();

  // initialize counter for dataset passes
  size_t cntDataPasses = 0;

  // refinement variables
  size_t refNum = adaptivityConfig.numRefinements_;
  size_t numPoints = adaptivityConfig.noPoints_;
  size_t refCnt = 0;
  double threshold = adaptivityConfig.threshold_;
  double currentBatchError = 0.0;
  double currentTrainError = 0.0;
  bool doRefine = false;  // set true by convergence monitor
  std::shared_ptr<ConvergenceMonitor> monitor(new ConvergenceMonitor(
      errorDeclineThreshold, errorDeclineBufferSize, minRefInterval));

  // decay factor for addaptive learning rate
  double lGamma = 0.001;

  // auxiliary variable for accuracy (error) measurement
  double acc = getAccuracy(*testData, *testLabels, 0.0);
  avgErrors.append(1.0 - acc);

  // parameters for ADAM
  /*sgpp::base::DataVector m(alpha->getSize(), 0.0);
  sgpp::base::DataVector v(alpha->getSize(), 0.0);
  double beta_1 = 0.25;
  double beta_2 = 0.95;
  double epsilon = 1e-8;*/

  // counts total number of processed data points
  size_t processedPoints = 0;
  // main loop which performs the learning process
  while (cntDataPasses < maxDataPasses) {
    for (size_t currIt = 0; currIt < trainData->getNrows(); currIt++) {
      // get next training sample x and its label y
      sgpp::base::DataVector x(dim);
      trainData->getRow(currIt, x);
      double y = trainLabels->get(currIt);

      // store data point in batch dataset used for checking
      // predictive refinement criterion
      // if validation set is used -> not needed
      if (!useValidData) {
        pushToBatch(x, y);
      }

      sgpp::base::DataVector delta(alpha->getSize());

      sgpp::base::DataVector singleAlpha(1);
      singleAlpha[0] = 1.0;

      // perform SGD step
      sgpp::base::DataMatrix dm(x.getPointer(), 1, x.getSize());
      std::unique_ptr<base::OperationMultipleEval> multEval(
          op_factory::createOperationMultipleEval(*grid, dm));
      multEval->multTranspose(singleAlpha, delta);

      double residual = delta.dotProduct(*alpha) - y;

      // ADAM
      // gradient
      /*delta.mult(residual);
      delta.axpy(lambda, *alpha);
      sgpp::base::DataVector grad_1(delta);
      sgpp::base::DataVector grad_2(delta);
      //update biased first and second moment estimates
      grad_1.mult(1.0-beta_1);
      m.mult(beta_1);
      m.add(grad_1);
      grad_2.sqr();
      grad_2.mult(1.0-beta_2);
      v.mult(beta_2);
      v.add(grad_2);
      //update bias-corrected first and second moment estimates
      //m.mult(1.0/(1.0-std::pow(beta_1,currIt+1)));
      //v.mult(1.0/(1.0-std::pow(beta_2,currIt+1)));
      m.mult(1.0/(1.0-std::pow(beta_1,processedPoints+1)));
      v.mult(1.0/(1.0-std::pow(beta_2,processedPoints+1)));
      //update alpha
      v.sqrt();
      sgpp::base::DataVector eps_vec(alpha->getSize(), epsilon);
      v.add(eps_vec);
      m.componentwise_div(v);
      m.mult(currentGamma);
      alpha->sub(m);*/

      // SGD
      alpha->mult(1 - currentGamma * lambda);
      alpha->axpy(-currentGamma * residual, delta);

      // learning rate according to L. Bottou
      currentGamma =
          gamma *
          std::pow(
              (1 + gamma * lGamma * (static_cast<double>(processedPoints) + 1)),
              -0.75);
      // currentGamma = gamma * std::pow(
      //  (1 + gamma*lambda*(static_cast<double>(processedPoints)+1)), -0.75);

      // smoothing according to L. Bottou
      size_t t1 = (processedPoints > dim + 1) ? processedPoints - dim : 1;
      size_t t2 = (processedPoints > trainData->getNrows() + 1)
                      ? processedPoints - trainData->getNrows()
                      : 1;
      double mu = (t1 > t2) ? static_cast<double>(t1) : static_cast<double>(t2);
      mu = 1.0 / mu;

      // average SGD / ADAM
      alphaAvg->mult(1 - mu);
      alphaAvg->axpy(mu, *alpha);

      // check if refinement should be performed
      if (refMonitor == "periodic") {
        // periodic monitor
        if ((refCnt < refNum) && (processedPoints > 0) &&
            ((processedPoints + 1) % refPeriod == 0)) {
          doRefine = true;
        }
      } else if (refMonitor == "convergence") {
        // check convergence monitor
        if (refCnt < refNum) {
          currentBatchError = getError(*batchData, *batchLabels, "MSE");
          currentTrainError = getError(*trainData, *trainLabels, "MSE");
          monitor->pushToBuffer(currentBatchError, currentTrainError);
          if (monitor->nextRefCnt > 0) {
            monitor->nextRefCnt--;
          }
          if (monitor->nextRefCnt == 0) {
            doRefine = monitor->checkConvergence();
          }
        }
      }
      if (doRefine) {
        // acc = getAccuracy(*testData, *testLabels, 0.0);
        // avgErrors.append(1.0 - acc);
        std::cout << "refinement at iteration: " << processedPoints
                  << std::endl;

        base::GridStorage& gridStorage = grid->getStorage();

        HashRefinement refinement;

        if (refType == "predictive") {
          // predictive refinement based on error contributions
          PredictiveRefinement decorator(&refinement);
          getBatchError(*batchData, *batchLabels, *batchError);
          PredictiveRefinementIndicator indicator(*grid, *batchData,
                                                  *batchError, numPoints);
          decorator.free_refine(gridStorage, indicator);
        } else if (refType == "impurity") {
          // impurity-based refinement
          ImpurityRefinement decorator(&refinement);
          sgpp::base::DataVector predictedLabels(batchData->getNrows());
          predict(*batchData, predictedLabels);
          ImpurityRefinementIndicator indicator(
              *grid, *batchData, nullptr, nullptr, nullptr, predictedLabels,
              threshold, numPoints);
          decorator.free_refine(gridStorage, indicator);
        }
        alpha->resizeZero(grid->getSize());
        alphaAvg->resizeZero(grid->getSize());

        // required for ADAM
        // m.resizeZero(grid->getSize());
        // v.resizeZero(grid->getSize());

        std::cout << "refinement step: " << refCnt + 1 << std::endl;
        std::cout << "new grid size: " << grid->getSize() << std::endl;

        refCnt++;
        doRefine = false;
        if (refMonitor == "convergence") {
          monitor->nextRefCnt = monitor->minRefInterval;
        }
      }

      // save current error
      if ((processedPoints + 1) % 10 == 0) {
        acc = getAccuracy(*testData, *testLabels, 0.0);
        avgErrors.append(1.0 - acc);
      }

      processedPoints++;
    }
    cntDataPasses++;
  }
  std::cout << "# Training finished" << std::endl;
  std::cout << "final grid size: " << grid->getSize() << std::endl;
  // double mse = getError(*testData, *testLabels, "MSE");
  // std::cout << "MSE: " << mse << std::endl;

  error = 1.0 - getAccuracy(*testData, *testLabels, 0.0);
}

void LearnerSGD::storeResults(base::DataMatrix& testDataset) {
  base::DataVector predictedLabels(testDataset.getNrows());
  predict(testDataset, predictedLabels);

  std::ofstream output;
  // write predicted class labels to csv file
  output.open("ASGD_predicted_classes.csv");
  if (output.fail()) {
    std::cout << "failed to create csv file!" << std::endl;
  } else {
    for (size_t i = 0; i < predictedLabels.getSize(); i++) {
      base::DataVector x(2);
      testDataset.getRow((size_t)i, x);
      output << x[0] << ";" << x[1] << ";" << predictedLabels[i] << std::endl;
    }
    output.close();
  }
  // write grid to csv file
  output.open("ASGD_grid.csv");
  if (output.fail()) {
    std::cout << "failed to create csv file!" << std::endl;
  } else {
    base::GridStorage& storage = grid->getStorage();
    base::GridStorage::grid_map_iterator end_iter = storage.end();
    for (base::GridStorage::grid_map_iterator iter = storage.begin();
         iter != end_iter; iter++) {
      base::DataVector gpCoord(testDataset.getNcols());
      storage.getCoordinates(*(iter->first), gpCoord);
      for (size_t d = 0; d < gpCoord.getSize(); d++) {
        if (d < gpCoord.getSize() - 1) {
          output << gpCoord[d] << ";";
        } else {
          output << gpCoord[d] << std::endl;
        }
      }
    }
    output.close();
  }

  // write function evaluations to csv file
  double stepSize = 0.01;
  base::DataMatrix values(0, 2);
  base::DataVector range(101);
  for (size_t i = 0; i < 101; i++) {
    range.set(i, stepSize * (static_cast<double>(i)));
  }
  for (size_t i = 0; i < range.getSize(); i++) {
    for (size_t j = 0; j < range.getSize(); j++) {
      base::DataVector row(2);
      row.set(1, range.get(i));
      row.set(0, range.get(j));
      values.appendRow(row);
    }
  }
  // evaluate learned function at all points from values
  // and write result to csv file
  output.open("ASGD_fun_evals.csv");
  std::unique_ptr<base::OperationEval> opEval(
      op_factory::createOperationEval(*grid));
  for (size_t i = 0; i < values.getNrows(); i++) {
    // get next test sample x
    base::DataVector x(2);
    values.getRow(i, x);
    double res = opEval->eval(*alphaAvg, x);
    output << res << ";" << std::endl;
  }
  output.close();
}

double LearnerSGD::getError(sgpp::base::DataMatrix& data,
                            sgpp::base::DataVector& labels,
                            std::string errorType) {
  size_t numData = data.getNrows();
  sgpp::base::DataVector result(numData);
  sgpp::base::DataVector error(numData);
  error.setAll(0.0);

  std::unique_ptr<base::OperationMultipleEval> opEval(
      op_factory::createOperationMultipleEval(*grid, data));
  opEval->mult(*alphaAvg, result);

  double res = -1.0;
  if (errorType == "MSE") {
    for (size_t i = 0; i < numData; i++) {
      error.set(i, labels.get(i) - result.get(i));
    }
    // MSE
    double sum = 0;
    for (size_t i = 0; i < numData; i++) {
      sum += error.get(i) * error.get(i);
    }
    res = (sum / static_cast<double>(numData));
  }
  if (errorType == "Hinge") {
    for (size_t i = 0; i < numData; i++) {
      error.set(i, std::max(0.0, 1.0 - labels.get(i) * result.get(i)));
    }
    // Hinge
    double sum = 0;
    for (size_t i = 0; i < numData; i++) {
      sum += error.get(i);
    }
    res = (sum / static_cast<double>(numData));
  }
  return res;
}

void LearnerSGD::getBatchError(sgpp::base::DataMatrix& data,
                               sgpp::base::DataVector& labels,
                               sgpp::base::DataVector& error) {
  size_t numData = data.getNrows();
  sgpp::base::DataVector result(numData);

  std::unique_ptr<base::OperationMultipleEval> opEval(
      op_factory::createOperationMultipleEval(*grid, data));
  opEval->mult(*alphaAvg, result);

  for (size_t i = 0; i < numData; i++) {
    error.set(i, pow(labels.get(i) - result.get(i), 2));
  }
}

double LearnerSGD::getAccuracy(sgpp::base::DataMatrix& testData,
                               sgpp::base::DataVector& testLabels,
                               double threshold) {
  sgpp::base::DataVector predictedLabels(testData.getNrows());
  predict(testData, predictedLabels);

  return getAccuracy(testLabels, threshold, predictedLabels);
}

double LearnerSGD::getAccuracy(sgpp::base::DataVector& testLabels,
                               double threshold,
                               sgpp::base::DataVector& predictedLabels) {
  double result = -1.0;

  if (predictedLabels.getSize() != testLabels.getSize()) {
    throw base::application_exception(
        "LearnerSGD::getAccuracy: lengths of label-vectors do not match!");
  }

  size_t correct = 0;
  for (size_t i = 0; i < predictedLabels.getSize(); i++) {
    if ((predictedLabels.get(i) >= threshold && testLabels.get(i) >= 0.0) ||
        (predictedLabels.get(i) < threshold && testLabels.get(i) < 0.0)) {
      correct++;
    }
  }

  result = static_cast<double>(correct) /
           static_cast<double>(predictedLabels.getSize());

  return result;
}

void LearnerSGD::predict(sgpp::base::DataMatrix& testData,
                         sgpp::base::DataVector& predictedLabels) {
  predictedLabels.resize(testData.getNrows());
  sgpp::base::DataVector result(testData.getNrows());

  std::unique_ptr<base::OperationMultipleEval> opEval(
      op_factory::createOperationMultipleEval(*grid, testData));
  opEval->mult(*alphaAvg, result);

  for (size_t i = 0; i < testData.getNrows(); i++) {
    if (result.get(i) >= 0.0) {
      predictedLabels.set(i, 1.0);
    } else {
      predictedLabels.set(i, -1.0);
    }
  }
}

void LearnerSGD::pushToBatch(sgpp::base::DataVector& x, double y) {
  static size_t nextIdx = 0;
  if (batchData->getUnused() > 0) {
    batchData->appendRow(x);
    (*batchLabels)[nextIdx] = y;
  } else {
    batchData->setRow(nextIdx, x);
    (*batchLabels)[nextIdx] = y;
  }
  nextIdx = (nextIdx + 1) % batchSize;
}

}  // namespace datadriven
}  // namespace sgpp

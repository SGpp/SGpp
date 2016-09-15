// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/LearnerSGD.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/PredictiveRefinement.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/functors/PredictiveRefinementIndicator.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>
#include <cmath>

using sgpp::base::GridStorage;
using sgpp::base::HashRefinement;
using sgpp::base::PredictiveRefinement;
using sgpp::base::PredictiveRefinementIndicator;

namespace sgpp {
namespace datadriven {

LearnerSGD::LearnerSGD(sgpp::base::RegularGridConfiguration& gridConfig,
                       sgpp::base::AdpativityConfiguration& adaptivityConfig
                       /*sgpp::datadriven::RegularizationType& regularization*/) 
    : gridConfig(gridConfig),
      adaptivityConfig(adaptivityConfig),
      grid(nullptr),
      alpha(nullptr),
      alphaAvg(nullptr),
      trainData(nullptr),
      trainLabels(nullptr),
      testData(nullptr),
      testLabels(nullptr),
      batchData(nullptr),
      batchLabels(nullptr),
      batchError(nullptr),
      currentBatchError(0),
      lambda(0),
      gamma(0),
      currentGamma(0),
      batchSize(0),
      smoothedErrorDecline(0),
      smoothedErrorDeclineBufferSize(0),
      error(0) {}


LearnerSGD::~LearnerSGD() {}


void LearnerSGD::initialize(sgpp::base::DataMatrix& pTrainData,
                            sgpp::base::DataVector& pTrainLabels,
                            sgpp::base::DataMatrix& pTestData,
                            sgpp::base::DataVector& pTestLabels,
                            double lambda,
                            double gamma,
                            double smoothedErrorDecline,
                            size_t batchSize,
                            size_t bufferSize) {
  trainData = std::make_shared<base::DataMatrix>(pTrainData);
  trainLabels = std::make_shared<base::DataVector>(pTrainLabels);
  testData = std::make_shared<base::DataMatrix>(pTestData);
  testLabels = std::make_shared<base::DataVector>(pTestLabels);

  this->lambda = lambda;
  this->gamma = gamma;
  currentGamma = gamma;

  this->batchSize = batchSize;

  //batchData = new base::DataMatrix(0, trainData->getNcols());
  batchData = std::shared_ptr<base::DataMatrix>(new base::DataMatrix(0, trainData->getNcols()));
  batchData->addSize(batchSize);
  batchData->setAll(0.0);
  //batchLabels = new base::DataVector(batchSize);
  batchLabels = std::shared_ptr<base::DataVector>(new base::DataVector(batchSize));
  batchLabels->setAll(0.0);
  //batchError = new base::DataVector(batchSize);
  batchError = std::shared_ptr<base::DataVector>(new base::DataVector(batchSize));
  batchError->setAll(0.0);

  smoothedErrorDeclineBufferSize = bufferSize;
  this->smoothedErrorDecline = smoothedErrorDecline;

  gridConfig.dim_ = trainData->getNcols();
  grid = createRegularGrid();
  std::cout << "# grid points: " << grid->getSize() << std::endl;

  //alpha = new base::DataVector(grid->getSize());
  alpha = std::shared_ptr<base::DataVector>(new base::DataVector(grid->getSize()));
  alpha->setAll(0.0);
  //alphaAvg = new base::DataVector(grid->getSize());
  alphaAvg = std::shared_ptr<base::DataVector>(new base::DataVector(grid->getSize()));
  alphaAvg->setAll(0.0);
}


std::shared_ptr<base::Grid> LearnerSGD::createRegularGrid() {
  // load grid
  std::unique_ptr<base::Grid> uGrid;
  if (gridConfig.type_ == base::GridType::Linear) {
    uGrid = base::Grid::createLinearGrid(gridConfig.dim_);
  } 
  else if (gridConfig.type_ == base::GridType::ModLinear) {
    uGrid = base::Grid::createModLinearGrid(gridConfig.dim_);
  } 
  else {
    throw base::application_exception("LearnerSGD::initialize : grid type is not supported");
  }

  uGrid->getGenerator().regular(gridConfig.level_);

  // move the grid to be shared
  std::shared_ptr<base::Grid> sGrid{std::move(uGrid)};

  return sGrid;
}


void LearnerSGD::train(size_t dataNum) {
  size_t dim = trainData->getNcols();  
  size_t maxRuns = 2; //ToDo: pass to learner as parameter
  size_t numRuns = 0;

  size_t refNum = adaptivityConfig.numRefinements_;
  size_t numPoints = adaptivityConfig.noPoints_;
  size_t refSteps = 0;
  double threshold = adaptivityConfig.threshold_;
  size_t refPeriod = 15;

  // refinement monitor variables
  double oldErrorSum = 0.;
  double oldErrorLast = 0.;
  currentBatchError = 0.;
  double ratio = 1.;

  // for error plotting
  sgpp::base::DataVector errors;
  //sgpp::base::DataVector grids;
  // for plotting only
  //double acc = getAccuracy(*trainData, *trainLabels, 0.0);
  double acc = getAccuracy(*testData, *testLabels, 0.0);
  errors.append(1.0 - acc);

  while (numRuns < maxRuns) {			
    for (size_t curr_it = 0; curr_it < trainData->getNrows(); curr_it++) {
      // get next training sample x and its label y
      sgpp::base::DataVector x(dim);					
      trainData->getRow(curr_it, x);
      double y = trainLabels->get(curr_it);  
      
      // store in batch dataset used for checking convergence
      pushToBatch(x, y);
      
      sgpp::base::DataVector delta(alpha->getSize());

      sgpp::base::DataVector singleAlpha(1);
      singleAlpha[0] = 1.0;

      // perform SGD step
      sgpp::base::DataMatrix dm(x.getPointer(), 1, x.getSize());
      //sgpp::base::OperationMultipleEval* multEval = 
      //  sgpp::op_factory::createOperationMultipleEval(*grid, dm);
      std::unique_ptr<base::OperationMultipleEval> multEval(op_factory::createOperationMultipleEval(*grid, dm));
      multEval->multTranspose(singleAlpha, delta);
      //delete multEval;

      double residual = delta.dotProduct(*alpha) - y;

      alpha->mult(1-currentGamma * lambda);
      alpha->axpy(-currentGamma * residual, delta);

      //double mu = 0.1; // boring exponential smoothing

      // L. Bottou exciting smoothing
      size_t t1 = (curr_it > dim+1) ? curr_it - dim : 1;
      size_t t2 = (curr_it > trainData->getNrows()+1) ? curr_it - trainData->getNrows() : 1;
      double mu = (t1>t2) ? static_cast<double>(t1) : static_cast<double>(t2);
      mu = 1.0/mu;

      alphaAvg->mult(1-mu);
      alphaAvg->axpy(mu, *alpha);

      // Calculate smoothed error
      //currentBatchError = getError(*batchData, *batchLabels, "MSE");
      if (refSteps < refNum) {
        currentBatchError = getError(*testData, *testLabels, "MSE");
      }
      if (smoothedErrorDeclineBuffer.size() >= smoothedErrorDeclineBufferSize) {
        // Calculate average of old minibatch errors
        for (std::list<double>::iterator it = smoothedErrorDeclineBuffer.begin();
             it != smoothedErrorDeclineBuffer.end();
             ++it) {
          oldErrorSum += *it;
        }

        oldErrorLast = smoothedErrorDeclineBuffer.back();

        // Update errorOnMinibatch
        smoothedErrorDeclineBuffer.pop_back();
        smoothedErrorDeclineBuffer.push_front(currentBatchError);

        // Update ratio
        ratio = (oldErrorLast - currentBatchError) / oldErrorSum;
        ratio *= 1.0 / (double)smoothedErrorDeclineBufferSize;

      } 
      else {
        smoothedErrorDeclineBuffer.push_front(currentBatchError);
      }

      // refinement
      if ( (refSteps < refNum) && (curr_it > 0 ) && ((curr_it+1) % refPeriod == 0) ) {
      //if ( (refSteps < refNum) && (ratio <= smoothedErrorDecline) ) {
        //std::cout << "threshold: " << smoothedErrorDecline << std::endl;
        //std::cout << "iteration: " << curr_it << std::endl;
        //std::cout << "ratio: " << ratio << std::endl;
        double acc_ref = getAccuracy(*testData, *testLabels, 0.0);
        std::cout << "accuracy before refinement: " << acc_ref << std::endl;
        // Predictive refinement based on error measurements
        base::GridStorage& gridStorage = grid->getStorage();
        HashRefinement refinement;
        PredictiveRefinement decorator(&refinement);
        
        getBatchError(*batchData, *batchLabels, *batchError);
        //std::cout << "Error over all = "  << batchError->sum() << std::endl;
        PredictiveRefinementIndicator indicator(*grid, *batchData, *batchError, numPoints);
        decorator.free_refine(gridStorage, indicator);

        alpha->resizeZero(grid->getSize());
        alphaAvg->resizeZero(grid->getSize());

        std::cout << "refinement step: " << refSteps+1 << std::endl;
        std::cout << "new grid size: " << grid->getSize() << std::endl;

        oldErrorSum = 0.;
        oldErrorLast = 0.;
        currentBatchError = 0.;
        ratio = 1.;

        refSteps++;      
      }
      acc = getAccuracy(*testData, *testLabels, 0.0);
      errors.append(1.0 - acc);
    }
    numRuns++;
  }
  double mse = getError(*testData, *testLabels, "MSE");
  std::cout << "MSE: " << mse << std::endl;
  acc = getAccuracy(*testData, *testLabels, 0.0);
  std::cout << "accuracy: " << acc << std::endl;

  //write error evaluation to .csv
  std::ofstream output;
  output.open("SGD_ripley_err_rate_predictive_train_"+std::to_string(dataNum)+".csv");
  //output.open("ripley_err_rate_measure_train_1.csv");
  //output.open("ripley_err_rate_impurity_train_1.csv");
  //output.open("banana_err_rate_measure_train_1.csv");
  //output.open("banana_err_rate_impurity_train_1.csv");
  if (output.fail()) {
    std::cout << "failed to create .csv file!" << std::endl;  
  }
  else {
    for (size_t i = 0; i < errors.getSize(); i++) {					
      output << errors.get(i) << ";" << std::endl;
    }
    output.close();
  }

  error = 1.0 - getAccuracy(*testData, *testLabels, 0.0);
}


void LearnerSGD::storeResults(base::DataMatrix& testDataset,
                              base::DataVector& testLabels,
                              double threshold) {

  base::DataVector computedLabels(testDataset.getNrows());
  predict(testDataset, computedLabels);

  std::ofstream output;
  //write computed labels to .csv
  //output.open("banana_SGD_predicted_train_1.csv");
  //output.open("banana_SGD_predicted_train_1.csv");
  output.open("ripley_SGD_predicted.csv");
  if (output.fail()) {
    std::cout << "failed to create .csv file!" << std::endl;  
  }
  else {
    for (size_t i = 0; i < computedLabels.getSize(); i++) {
      base::DataVector x(2);					
      testDataset.getRow((size_t)i, x);
      output << x[0] << ";" << x[1] << ";" << computedLabels[i] << std::endl;
    }
    output.close();
  }
  //write grid to .csv
  //output.open("banana_impurity_grid_train_1.csv");
  //output.open("banana_surplus_grid_train_1.csv");
  output.open("ripley_SGD_grid.csv");
  if (output.fail()) {
    std::cout << "failed to create .csv file!" << std::endl;  
  }
  else {
    base::GridStorage& storage = grid->getStorage();
    base::GridStorage::grid_map_iterator end_iter = storage.end();
    for (base::GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter; iter++) { 
      base::DataVector gpCoord(testDataset.getNcols());
      storage.getCoordinates(*(iter->first), gpCoord);
      for (size_t d = 0; d < gpCoord.getSize(); d++) {
        if (d < gpCoord.getSize()-1) {
          output << gpCoord[d] << ";";
        }
        else {
          output << gpCoord[d] << std::endl;
        }
      }
    }
    output.close();
  }

  //write function evaluations to .csv
  double stepSize = 0.01;
  base::DataMatrix values(0,2);
  //std::cout << values.getNrows() << std::endl;
  base::DataVector range(101);
  for (size_t i = 0; i < 101; i++) {
    range.set(i, stepSize*(static_cast<double>(i)));
  }
  //std::cout << range.getSize() << std::endl;
  for (size_t i = 0; i < range.getSize(); i++) {
    for (size_t j = 0; j < range.getSize(); j++) {
      base::DataVector row(2);
      row.set(1, range.get(i));
      row.set(0, range.get(j));
      values.appendRow(row);
    }
  }
  //std::cout << values.getNrows() << std::endl;
  // evaluate learned function at all points from values
  // and write result to csv file
  output.open("SGD_fun_evals.csv");
  std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEval(*grid));
  for (size_t i = 0; i < values.getNrows(); i++) {
    // Get next test sample x 
    base::DataVector x(2);					
    values.getRow(i, x);
    //std::cout << x[0] << " , " << x[1] << std::endl;
    double res = opEval->eval(*alphaAvg, x);
    output << res << ";" << std::endl;
  }  
  output.close();

}


double LearnerSGD::getError(sgpp::base::DataMatrix& data, sgpp::base::DataVector& labels, 
                            std::string errorType) {
  size_t numData = data.getNrows();
  sgpp::base::DataVector result(numData);
  sgpp::base::DataVector error(numData);
  error.setAll(0.0);

  //OperationMultipleEval* eval = sgpp::op_factory::createOperationMultipleEval(*grid, data);
  std::unique_ptr<base::OperationMultipleEval> opEval(op_factory::createOperationMultipleEval(*grid, data));
  opEval->mult(*alphaAvg, result);

  double res = -1.0;
  if (errorType == "MSE") {
    for (size_t i = 0; i < numData; i++) {
            error.set(i, labels.get(i) - result.get(i));
    }
    // Error
    double sum = 0;
    for (size_t i = 0; i < numData; i++) {
      sum += error.get(i) * error.get(i);
    }
    res = (sum / (double)numData);
  }
  /*if (errorType == "ACCURACY") {
    size_t correct = 0;
    for (size_t i = 0; i < numData; i++) {
      correct += (result.get(i) < 0) == (labels.get(i) < 0) ? 1 : 0;
      error.set(i, labels.get(i) - result.get(i));
    }
    res = static_cast<double>(correct) / static_cast<double>(numData);
  }*/
  return res;
}


void LearnerSGD::getBatchError(sgpp::base::DataMatrix& data, sgpp::base::DataVector& labels, 
                               sgpp::base::DataVector& error) {
  size_t numData = data.getNrows();
  sgpp::base::DataVector result(numData);

  //OperationMultipleEval* eval = sgpp::op_factory::createOperationMultipleEval(*grid, data);
  std::unique_ptr<base::OperationMultipleEval> opEval(op_factory::createOperationMultipleEval(*grid, data));
  opEval->mult(*alphaAvg, result);
 
  for (size_t i = 0; i < numData; i++) {
    error.set(i, pow(labels.get(i) - result.get(i), 2));
  }
}


double LearnerSGD::getAccuracy(sgpp::base::DataMatrix& testData, sgpp::base::DataVector& testLabels, 
                               double threshold) {
  sgpp::base::DataVector computedLabels(testData.getNrows());
  predict(testData, computedLabels);

  return getAccuracy(computedLabels, testLabels, threshold);
}


double LearnerSGD::getAccuracy(sgpp::base::DataVector& computedLabels,
                               sgpp::base::DataVector& testLabels,
                               double threshold) {
  double result = -1.0;

  if (computedLabels.getSize() != testLabels.getSize()) {
    throw base::application_exception(
        "LearnerSGD::getAccuracy: lengths of classes vectors do not match!");
  }

  size_t correct = 0;
  for (size_t i = 0; i < computedLabels.getSize(); i++) {
    if ((computedLabels.get(i) >= threshold && testLabels.get(i) >= 0.0) ||
       (computedLabels.get(i) < threshold && testLabels.get(i) < 0.0)) {
      correct++;
    }
  }

  result = static_cast<double>(correct) / static_cast<double>(computedLabels.getSize());

  return result;
}

void LearnerSGD::predict(sgpp::base::DataMatrix& testData,
                         sgpp::base::DataVector& computedLabels) {
  computedLabels.resize(testData.getNrows());
  sgpp::base::DataVector result(testData.getNrows());
  size_t dim = testData.getNcols();  

  std::unique_ptr<base::OperationMultipleEval> opEval(op_factory::createOperationMultipleEval(*grid, testData));
  opEval->mult(*alphaAvg, result);

  for (size_t i = 0; i < testData.getNrows(); i++) {
    if (result.get(i) >= 0.0) {
      computedLabels.set(i, 1.0); 
    }
    else {
      computedLabels.set(i, -1.0);
    }
  }
}


void LearnerSGD::pushToBatch(sgpp::base::DataVector& x, double y) {
  static size_t next_idx = 0;
  if (batchData->getUnused() > 0) {
    batchData->appendRow(x);
    (*batchLabels)[next_idx] = y;
  }
  else {
    batchData->setRow(next_idx, x);
    (*batchLabels)[next_idx] = y;
  }
  next_idx = (next_idx + 1) % batchSize;
}


/*int LearnerSGD::getRandom(int limit) {
  int divisor = RAND_MAX / (limit + 1);
  int r;

  do {
    r = rand() / divisor;
  } while (r > limit);

  return r;
}*/


}  // namespace datadriven
}  // namespace sgpp


// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/LearnerSGD.hpp>
#include <sgpp/datadriven/algorithm/ConvergenceMonitor.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/PredictiveRefinement.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/functors/PredictiveRefinementIndicator.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/ImpurityRefinement.hpp>
#include <sgpp/base/grid/generation/functors/ImpurityRefinementIndicator.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>
#include <cmath>

using sgpp::base::GridStorage;
using sgpp::base::HashRefinement;
using sgpp::base::PredictiveRefinement;
using sgpp::base::PredictiveRefinementIndicator;
using sgpp::base::ImpurityRefinement;
using sgpp::base::ImpurityRefinementIndicator;

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
                            std::shared_ptr<sgpp::base::DataMatrix> pValData,
                            std::shared_ptr<sgpp::base::DataVector> pValLabels,
                            double lambda,
                            double gamma,
                            double smoothedErrorDecline,
                            size_t batchSize,
                            size_t bufferSize,
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

  // for periodic refinement -> don't use validation set for error measurement
  if (not useValidData) {
    batchData = std::shared_ptr<base::DataMatrix>(new base::DataMatrix(0, trainData->getNcols()));
    batchData->addSize(batchSize);
    batchData->setAll(0.0);
    batchLabels = std::shared_ptr<base::DataVector>(new base::DataVector(batchSize));
    batchLabels->setAll(0.0);
  }
  else {
    batchData = pValData;
    batchLabels = pValLabels;
  }

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
  size_t maxRuns = 15; //ToDo: pass to learner as parameter
  size_t numRuns = 0;

  size_t refNum = adaptivityConfig.numRefinements_;
  size_t numPoints = adaptivityConfig.noPoints_;
  size_t refSteps = 0;
  double threshold = adaptivityConfig.threshold_;

  // refinement monitor variables
  // if periodic refinement is used
  size_t refPeriod = 200;
  // if error based refinement is used
  double currentBatchError = 0.0;
  double currentTrainError = 0.0;
  bool doRefine = false;
  std::shared_ptr<ConvergenceMonitor> monitor(new ConvergenceMonitor(smoothedErrorDecline,
                                                                     smoothedErrorDeclineBufferSize,
                                                                     100));

  double mse;

  // decay factor for addaptive learning rate
  double l_gamma = 0.001;

  // for error plotting
  sgpp::base::DataVector errors;
  //sgpp::base::DataVector grids;
  // for plotting only
  double acc = getAccuracy(*testData, *testLabels, 0.0);
  //double acc = getAccuracy(*trainData, *trainLabels, 0.0);
  //errors.append(1.0 - acc);  //simple
  avgErrors.append(1.0 - acc);  //5-fold
  sgpp::base::DataVector runtimes;

  sgpp::base::DataVector m(alpha->getSize(), 0.0);
  sgpp::base::DataVector v(alpha->getSize(), 0.0);
  double beta_1 = 0.005; // 0.9
  double beta_2 = 0.95; // 0.999
  double epsilon = 1e-8; //1e-8

  size_t processedPoints = 0;
  while (numRuns < maxRuns) {		
    ////////////////////////////////////////////////
    /*    bool bufferFull;
        base::DataMatrix* validBatch = new base::DataMatrix(0,4);
        base::DataMatrix* trainBatch = new base::DataMatrix(0,4);
        base::DataVector* validBatch_c = new base::DataVector(0);
        base::DataVector* trainBatch_c = new base::DataVector(0);
        for (size_t idx_ = 0; idx_ < 250; idx_++) { //250
          int divisor = RAND_MAX / static_cast<int>(batchData->getNrows());
          int r;
          do {
            r = rand() / divisor;
          } while (r > (static_cast<int>(batchData->getNrows())-1));
          base::DataVector x(4);
          batchData->getRow((size_t)r, x);
          double y = batchLabels->get((size_t)r);
          validBatch->appendRow(x);
          validBatch_c->append(y);
        }
        for (size_t idx_ = 0; idx_ < 250; idx_++) {
          int divisor = RAND_MAX / static_cast<int>(trainData->getNrows());
          int r;
          do {
            r = rand() / divisor;
          } while (r > (static_cast<int>(trainData->getNrows())-1));
          base::DataVector x(4);
          trainData->getRow((size_t)r, x);
          double y = trainLabels->get((size_t)r);
          trainBatch->appendRow(x);
          trainBatch_c->append(y);
        }*/
    ///////////////////////////////////////////////////////////////////
    for (size_t curr_it = 0; curr_it < trainData->getNrows(); curr_it++) {
      // get next training sample x and its label y
      sgpp::base::DataVector x(dim);					
      trainData->getRow(curr_it, x);
      double y = trainLabels->get(curr_it);  
      
      // store datapoints in batch dataset used for checking predictive refinement criterion
      // if validation set is used -> not needed
      if (not useValidData) {
        pushToBatch(x, y);
      }

      sgpp::base::DataVector delta(alpha->getSize());

      sgpp::base::DataVector singleAlpha(1);
      singleAlpha[0] = 1.0;

      //time measure start
      clock_t begin;
      //if ( (processedPoints == 0) or (processedPoints % 1000 == 0) ) {
      begin = clock();
      //} 

      // perform SGD step
      sgpp::base::DataMatrix dm(x.getPointer(), 1, x.getSize());
      std::unique_ptr<base::OperationMultipleEval> multEval(op_factory::createOperationMultipleEval(*grid, dm));
      multEval->multTranspose(singleAlpha, delta);

      double residual = delta.dotProduct(*alpha) - y;

      /////ADAM/////
      //gradient
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
      //m.mult(1.0/(1.0-std::pow(beta_1,curr_it+1)));
      //v.mult(1.0/(1.0-std::pow(beta_2,curr_it+1)));
      m.mult(1.0/(1.0-std::pow(beta_1,processedPoints+1)));
      v.mult(1.0/(1.0-std::pow(beta_2,processedPoints+1)));
      //update alpha
      v.sqrt();
      sgpp::base::DataVector eps_vec(alpha->getSize(), epsilon);
      v.add(eps_vec);
      m.componentwise_div(v);
      m.mult(currentGamma);
      alpha->sub(m);*/
      //////////////

      /////SGD/////
      alpha->mult(1-currentGamma * lambda);
      alpha->axpy(-currentGamma * residual, delta);

      //std::cout << currentGamma << std::endl;
      //currentGamma = gamma * std::pow((1 + gamma*l_gamma*(static_cast<double>(curr_it)+1)), -0.75);
      currentGamma = gamma * std::pow((1 + gamma*l_gamma*(static_cast<double>(processedPoints)+1)), -0.75);

      //double mu = 0.1; // boring exponential smoothing

      // L. Bottou exciting smoothing
      //size_t t1 = (curr_it > dim+1) ? curr_it - dim : 1;
      //size_t t2 = (curr_it > trainData->getNrows()+1) ? curr_it - trainData->getNrows() : 1;
      size_t t1 = (processedPoints > dim+1) ? processedPoints - dim : 1;
      size_t t2 = (processedPoints > trainData->getNrows()+1) ? processedPoints - trainData->getNrows() : 1;
      double mu = (t1>t2) ? static_cast<double>(t1) : static_cast<double>(t2);
      mu = 1.0/mu;

      //ASGD - AADAM
      alphaAvg->mult(1-mu);
      alphaAvg->axpy(mu, *alpha);

      //time measure end
      /*clock_t end;
      if ( (processedPoints == 0) or (processedPoints % 1000 == 0) ) {
        end = clock();
        double elapsed_secs = double(end-begin)/CLOCKS_PER_SEC;
        //std::cout << "#time needed: " << elapsed_secs << std::endl;
      }*/

      //SGD
      //alphaAvg = alpha;

      //periodic monitor
      /*if ( (refSteps < refNum) && (processedPoints > 0 ) && ((processedPoints+1) % refPeriod == 0) ) {
        doRefine = true;
      }*/

      //convergence monitor
      if (refSteps < refNum) {
        currentBatchError = getError(*batchData, *batchLabels, "MSE");
        //currentBatchError = getError(*batchData, *batchLabels, "Hinge");
        //currentTrainError = getError(*trainData, *trainLabels, "MSE");
        //currentTrainError = getError(*trainData, *trainLabels, "Hinge");
        if (monitor->nextRefCnt > 0) {
          monitor->nextRefCnt--;
        }
        monitor->pushToBuffer(currentBatchError,currentTrainError);
        if ( (refSteps < refNum) && (monitor->nextRefCnt == 0) ){
          doRefine = monitor->checkConvergence();
        }
      }

      //convergence monitor
      /*if (refSteps < refNum) {
        /////////////////////////////////////////////////
        if (processedPoints % 500 == 0) { //%500
          bufferFull = false;
        }
        if (not bufferFull) {
        /////////////////////////////////////////////////
          //currentBatchError = getError(*validData, *validLabels, "MSE");
          //currentTrainError = getError(*trainData, *trainLabels, "MSE");
          currentBatchError = getError(*validBatch, *validBatch_c, "MSE");
          currentTrainError = getError(*trainBatch, *trainBatch_c, "MSE");
          monitor->pushToBuffer(currentBatchError,currentTrainError);
        ////
        }
        ////
        if (monitor->nextRefCnt > 0) {
          monitor->nextRefCnt--;
        }
        ///////////////////////////////////////////////
        if (monitor->validErrorDeclineBuffer.size() == smoothedErrorDeclineBufferSize) {
          bufferFull = true;
        /////////////////////////////////////////////////
          if (monitor->nextRefCnt == 0) {
            doRefine = monitor->checkConvergence();
            /////////////////////////////////////
            monitor->validErrorDeclineBuffer.clear();
            monitor->trainErrorDeclineBuffer.clear();
            /////////////////////////////////////
          }
        }
      }*/

      if (doRefine) {
        //acc = getAccuracy(*testData, *testLabels, 0.0);
        //avgErrors.append(1.0 - acc);
        std::cout << "refinement at iteration: " << processedPoints << std::endl;
        //std::cout << "ratio: " << batchRatio << std::endl;
        //std::cout << "mse (train data): " << getError(*trainData, *trainLabels, "MSE") << std::endl;

        // Predictive refinement based on error measurements
        base::GridStorage& gridStorage = grid->getStorage();

        HashRefinement refinement;

        PredictiveRefinement decorator(&refinement);
        getBatchError(*batchData, *batchLabels, *batchError);  
        PredictiveRefinementIndicator indicator(*grid, *batchData, *batchError, numPoints);

        /*ImpurityRefinement decorator(&refinement);
        sgpp::base::DataVector classesComputed(batchData->getNrows());
        predict(*batchData, classesComputed);
        ImpurityRefinementIndicator indicator(*grid, *batchData, nullptr, nullptr, nullptr, 
                                              classesComputed, threshold, numPoints);*/

        decorator.free_refine(gridStorage, indicator);

        alpha->resizeZero(grid->getSize());
        alphaAvg->resizeZero(grid->getSize());

        //required for ADAM
        //m.resizeZero(grid->getSize());	
        //v.resizeZero(grid->getSize());

        std::cout << "refinement step: " << refSteps+1 << std::endl;
        std::cout << "new grid size: " << grid->getSize() << std::endl;

        refSteps++;   
        doRefine = false;   

        monitor->nextRefCnt = monitor->minRefInterval;
      }

      clock_t end;
      //if ( (processedPoints == 0) or (processedPoints % 1000 == 0) ) {
      end = clock();
      double elapsed_secs = double(end-begin)/CLOCKS_PER_SEC;
      runtimes.append(elapsed_secs);
        //std::cout << "#time needed: " << elapsed_secs << std::endl;
      //}

      if ((processedPoints+1) % 50 == 0) {
      //if ((curr_it > 0) && ((curr_it+1) % 50 == 0)) {
        acc = getAccuracy(*testData, *testLabels, 0.0);
        avgErrors.append(1.0 - acc);  //5-fold
        //errors.append(1.0 - acc);  //simple
      }  
      /*acc = getAccuracy(*testData, *testLabels, 0.0);
      //avgErrors.append(1.0 - acc);  //5-fold
      errors.append(1.0 - acc);  //simple*/

      processedPoints++; 
    }
    std::cout << "final grid size: " << grid->getSize() << std::endl;

    numRuns++;
  }
  mse = getError(*testData, *testLabels, "MSE");
  std::cout << "MSE: " << mse << std::endl;
  /*acc = getAccuracy(*testData, *testLabels, 0.0);
  std::cout << "accuracy: " << acc << std::endl;*/

  //write error evaluation to .csv
  /*std::ofstream output;
  //output.open("ASGD_banana_error_monitor_predictive_runtimes_"+std::to_string(dataNum)+".csv");
  output.open("ASGD_banana_error_monitor_impurity_runtimes_"+std::to_string(dataNum)+".csv");
  if (output.fail()) {
    std::cout << "failed to create .csv file!" << std::endl;  
  }
  else {
    for (size_t i = 0; i < runtimes.getSize(); i++) {					
      output << runtimes.get(i) << ";" << std::endl;
    }
    output.close();
  }*/

  error = 1.0 - getAccuracy(*testData, *testLabels, 0.0);
  //error = 1.0 - getAccuracy(*batchData, *batchLabels, 0.0);
  //error = 1.0 - getAccuracy(*trainData, *trainLabels, 0.0);
  //error = getError(*testData, *testLabels, "MSE");
  //error = getError(*batchData, *batchLabels, "MSE");
  //error = getError(*trainData, *trainLabels, "MSE");
}


void LearnerSGD::storeResults(base::DataMatrix& testDataset,
                              base::DataVector& testLabels,
                              double threshold) {

  base::DataVector computedLabels(testDataset.getNrows());
  predict(testDataset, computedLabels);

  std::ofstream output;
  //write computed labels to .csv
  //output.open("ASGD_ripley_errormonitor_predictive_predicted.csv");
  //output.open("ASGD_ripley_periodicmonitor_impurity_predicted.csv");
  output.open("ASGD_banana_errormonitor_predictive_predicted.csv");
  //output.open("ASGD_banana_periodicmonitor_impurity_predicted.csv");
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
  //output.open("ASGD_ripley_errormonitor_predictive_grid.csv");
  //output.open("ASGD_ripley_periodicmonitor_impurity_grid.csv");
  output.open("ASGD_banana_errormonitor_predictive_grid.csv");
  //output.open("ASGD_banana_periodicmonitor_impurity_grid.csv");
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
  //output.open("ASGD_ripley_errormonitor_predictive_fun_evals.csv");
  //output.open("ASGD_ripley_periodicmonitor_impurity_fun_evals.csv");
  output.open("ASGD_banana_errormonitor_predictive_fun_evals.csv");
  //output.open("ASGD_banana_periodicmonitor_impurity_fun_evals.csv");
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
  if (errorType == "Hinge") {
    for (size_t i = 0; i < numData; i++) {
      error.set(i, std::max(0.0, 1.0 - labels.get(i) * result.get(i)));
    }
    // Error
    double sum = 0;
    for (size_t i = 0; i < numData; i++) {
      sum += error.get(i);
    }
    res = (sum / (double)numData);
  }
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


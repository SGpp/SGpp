// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/LearnerSVM.hpp>

#include <sgpp/datadriven/algorithm/ConvergenceMonitor.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/ForwardSelectorRefinement.hpp>
#include <sgpp/base/grid/generation/functors/ForwardSelectorRefinementIndicator.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/ImpurityRefinement.hpp>
#include <sgpp/base/grid/generation/functors/ImpurityRefinementIndicator.hpp>

#include <sgpp/base/exception/application_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

using sgpp::base::GridStorage;
using sgpp::base::HashRefinement;
using sgpp::base::ForwardSelectorRefinement;
using sgpp::base::ForwardSelectorRefinementIndicator;
using sgpp::base::ImpurityRefinement;
using sgpp::base::ImpurityRefinementIndicator;

namespace sgpp {
namespace datadriven {

// -------------------- constructors and destructors --------------------

LearnerSVM::LearnerSVM(sgpp::base::RegularGridConfiguration& gridConfig)                      
    : grid(nullptr),
      trainData(nullptr),
      trainLabels(nullptr),
      testData(nullptr),
      testLabels(nullptr),
      validData(nullptr),
      validLabels(nullptr),
      gridConfig(gridConfig),
      error(0.0) {}

LearnerSVM::~LearnerSVM() {}

// -----------------------------------------------------------------------

void LearnerSVM::initialize(base::DataMatrix& pTrainData, base::DataVector& pTrainLabels,
                            base::DataMatrix& pTestData, base::DataVector& pTestLabels,
                            std::shared_ptr<sgpp::base::DataMatrix> pValidData,
                            std::shared_ptr<sgpp::base::DataVector> pValidLabels) {
  trainData = std::make_shared<base::DataMatrix>(pTrainData);
  trainLabels = std::make_shared<base::DataVector>(pTrainLabels);
  testData = std::make_shared<base::DataMatrix>(pTestData);
  testLabels = std::make_shared<base::DataVector>(pTestLabels);
  validData = pValidData;
  validLabels = pValidLabels;
  gridConfig.dim_ = trainData->getNcols();
  grid = createRegularGrid();

  std::cout << "# grid points: " << grid->getSize() << std::endl;
  
  int B = 2400; //ToDo: pass budget B (=max number of stored support vectors) as parameter
               // DR10 16000
  
  // generate SVM
  svm = std::shared_ptr<PrimalDualSVM>(new PrimalDualSVM(grid->getSize(), trainData->getNcols(), B, false));
}

std::shared_ptr<base::Grid> LearnerSVM::createRegularGrid() {
  // load grid
  std::unique_ptr<base::Grid> uGrid;
  if (gridConfig.type_ == base::GridType::ModLinear) {
    uGrid = base::Grid::createModLinearGrid(gridConfig.dim_);
  } 
  else {
    throw base::application_exception("LearnerSVM::initialize : grid type is not supported");
  }

  uGrid->getGenerator().regular(gridConfig.level_);

  // move the grid to be shared
  std::shared_ptr<base::Grid> sGrid{std::move(uGrid)};

  return sGrid;
}

void LearnerSVM::train(size_t dataNum) {

  GridStorage& gridStorage = grid->getStorage();

  size_t dim = trainData->getNcols();
  
  size_t maxIterations = 15; //15
  size_t numIterations = 0;
  size_t t = 0;

  // refinement parameters - ToDo: pass to Learner as parameters (adaptivityConfig)
  size_t numRefSteps = 17; 
  size_t refSteps = 0; 
  size_t numPoints = 5; 
  double threshold = 0.0;
  double betaRef = 2.0;
 
  double lambda = 0.0015; 
  double eta = 0.0;
  double beta = 0.0;
  double rawPrediction = 0.0;
  
  //periodic monitor
  size_t refPeriod = 200;
  //convergence monitor
  double errorDeclineThreshold = 0.00005;
  size_t errorDeclineBufferSize = 300;
  double currentValidError = 0.0;
  double currentTrainError = 0.0;
  std::shared_ptr<ConvergenceMonitor> monitor(new ConvergenceMonitor(errorDeclineThreshold,
                                                                     errorDeclineBufferSize,
                                                                     100));
  bool doRefine = false;      // set true by monitor to trigger refinement

  // for error plotting
  sgpp::base::DataVector errors;
  //sgpp::base::DataVector grids;
  // for plotting only
  //double acc = getAccuracy(*trainData, *trainLabels, 0.0);
  double acc = getAccuracy(*testData, *testLabels, 0.0);
  //errors.append(1.0 - acc);  //simple
  avgErrors.append(1.0 - acc);  //5-fold
  sgpp::base::DataVector runtimes;

  size_t processedPoints = 0;
  size_t steps = 0;
    while (numIterations < maxIterations) {
      ////////////////////////////////////////////////
      /*  bool bufferFull;
        base::DataMatrix* validBatch = new base::DataMatrix(0,4);
        base::DataMatrix* trainBatch = new base::DataMatrix(0,4);
        base::DataVector* validBatch_c = new base::DataVector(0);
        base::DataVector* trainBatch_c = new base::DataVector(0);
        for (size_t idx_ = 0; idx_ < 250; idx_++) { //250
          int divisor = RAND_MAX / static_cast<int>(validData->getNrows());
          int r;
          do {
            r = rand() / divisor;
          } while (r > (static_cast<int>(validData->getNrows())-1));
          base::DataVector x(4);
          validData->getRow((size_t)r, x);
          double y = validLabels->get((size_t)r);
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
      for (size_t i = 0; i < trainData->getNrows(); i++) {
        // Get next training sample x and its label y
        //int idx = getRandom(static_cast<int>(trainData->getNrows()) - 1);
        sgpp::base::DataVector x(dim);					
        //trainData->getRow((size_t)idx, x);
        trainData->getRow(i, x);
        double y = trainLabels->get(i);   

        //time measure start
        clock_t begin;
        //if ( (i == 0) or (i % 1000 == 0) ) {
        begin = clock();
        //}                        
	  
        t += 1; 
        eta = 1.0 / (lambda * static_cast<double>(t));

        rawPrediction = svm->predictRaw(*grid, x, dim);
        svm->multiply(1.0 - lambda * eta);
        if (rawPrediction * y < 1.0) {
          beta = eta * y;
	  svm->add(*grid, x, beta, dim);
        }	

        //time measure end
        /*clock_t end;
        if ( (i == 0) or (i % 1000 == 0) ) {
          end = clock();
          double elapsed_secs = double(end-begin)/CLOCKS_PER_SEC;
          std::cout << "#time needed: " << elapsed_secs << std::endl;
        }*/

        //convergence monitor
        if (refSteps < numRefSteps) {
          currentValidError = getError(*validData, *validLabels, "Hinge");
          //currentValidError = getError(*validData, *validLabels, "Hinge");
          //currentTrainError = getError(*trainData, *trainLabels, "MSE");
          //currentTrainError = getError(*trainData, *trainLabels, "Hinge");
          if (monitor->nextRefCnt > 0) {
            monitor->nextRefCnt--;
          }
          monitor->pushToBuffer(currentValidError,currentTrainError);
          if ( (refSteps < numRefSteps) && (monitor->nextRefCnt == 0) ){
            doRefine = monitor->checkConvergence();
          }
        }

        //convergence monitor
      /*if (refSteps < numRefSteps) {
        /////////////////////////////////////////////////
        if (processedPoints % 500 == 0) { //%500
          bufferFull = false;
        }
        if (not bufferFull) {
        /////////////////////////////////////////////////
          //currentValidError = getError(*validData, *validLabels, "MSE");
          //currentTrainError = getError(*trainData, *trainLabels, "MSE");
          currentValidError = getError(*validBatch, *validBatch_c, "Hinge");
          currentTrainError = getError(*trainBatch, *trainBatch_c, "Hinge");
          monitor->pushToBuffer(currentValidError,currentTrainError);
        ////
        }
        ////
        if (monitor->nextRefCnt > 0) {
          monitor->nextRefCnt--;
        }
        ///////////////////////////////////////////////
        if (monitor->validErrorDeclineBuffer.size() == errorDeclineBufferSize) {
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
      
        //periodic
        /*if ( (refSteps < numRefSteps) && (steps > 0) && ((steps+1) % refPeriod == 0) ) {
          doRefine = true;
        }*/

        processedPoints ++;
        //std::cout << processedPoints << std::endl;
        if (doRefine) {
          //for plotting
          //acc = getAccuracy(*testData, *testLabels, 0.0);
          //avgErrors.append(1.0 - acc);
	  //
          std::cout << "refinement at iteration: " << steps+1 << std::endl;
          // compute current labels of support vectors 
          // (needed for impurity refinement)
          sgpp::base::DataVector svsClassesComputed((svm->svs)->getNrows());
          predict(*(svm->svs), svsClassesComputed);

          HashRefinement refinement;
          ForwardSelectorRefinement decorator(&refinement);
          //ImpurityRefinement decorator(&refinement);

          // refine points according to indicator
          ForwardSelectorRefinementIndicator indicator(*grid, *(svm->svs), *(svm->alphas), *(svm->w),
                                                       *(svm->w2), betaRef, threshold, numPoints);        
          //ImpurityRefinementIndicator indicator(*grid, *(svm->svs), (svm->alphas).get(), (svm->w).get(),
          //                                      (svm->w2).get(), svsClassesComputed, threshold, numPoints);

          decorator.free_refine(gridStorage, indicator);

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

        // for plotting only
        if (processedPoints % 50 == 0) { //banana %50 DR10 %1000
          acc = getAccuracy(*testData, *testLabels, 0.0);
          avgErrors.append(1.0 - acc);  //5-fold
          //errors.append(1.0 - acc);  //simple
          //output << grid->getSize() << ";" << std::endl;//new
        }
        //acc = getAccuracy(*testData, *testLabels, 0.0);
        //errors.append(1.0 - acc);  //simple
        //avgErrors.append(1.0 - acc);  //5-fold

        steps++;
      }
      std::cout << "final grid size: " << grid->getSize() << std::endl;
      double hinge = getError(*testData, *testLabels, "Hinge");
      std::cout << "hinge: " << hinge << std::endl;
      numIterations++;

    }
    //write error evaluation to .csv
    /*std::ofstream output;
    //output.open("SVM_banana_error_monitor_combinedmeasure_runtimes_"+std::to_string(dataNum)+".csv");
    output.open("SVM_banana_error_monitor_impurity_runtimes_"+std::to_string(dataNum)+".csv");
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

}

void LearnerSVM::storeResults(sgpp::base::DataMatrix& testDataset,
                              sgpp::base::DataVector& testLabels,
                              double threshold) {

  sgpp::base::DataVector computedLabels(testDataset.getNrows());
  predict(testDataset, computedLabels);

  std::ofstream output;
  //write computed labels to .csv
  //output.open("SVM_ripley_errormonitor_impurity_predicted.csv");
  //output.open("SVM_ripley_periodicmonitor_impurity_predicted.csv");
  output.open("SVM_banana_errormonitor_combinedmeasure_predicted.csv");
  //output.open("SVM_banana_periodicmonitor_combinedmeasure_predicted.csv");
  if (output.fail()) {
    std::cout << "failed to create .csv file!" << std::endl;  
  }
  else {
    for (size_t i = 0; i < computedLabels.getSize(); i++) {
      sgpp::base::DataVector x(2);					
      testDataset.getRow((size_t)i, x);
      output << x[0] << ";" << x[1] << ";" << computedLabels[i] << std::endl;
    }
    output.close();
  }
  //write grid to .csv
  //output.open("SVM_ripley_errormonitor_impurity_grid.csv");
  //output.open("SVM_ripley_periodicmonitor_impurity_grid.csv");
  output.open("SVM_banana_errormonitor_combinedmeasure_grid.csv");
  //output.open("SVM_banana_periodicmonitor_combinedmeasure_grid.csv");
  if (output.fail()) {
    std::cout << "failed to create .csv file!" << std::endl;  
  }
  else {
    sgpp::base::GridStorage& storage = grid->getStorage();
    sgpp::base::GridStorage::grid_map_iterator end_iter = storage.end();
    for (sgpp::base::GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter; iter++) { 
      sgpp::base::DataVector gpCoord(testDataset.getNcols());
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
  sgpp::base::DataMatrix values(0,2);
  sgpp::base::DataVector range(101);
  for (size_t i = 0; i < 101; i++) {
    range.set(i, stepSize*(static_cast<double>(i)));
  }
  for (size_t i = 0; i < range.getSize(); i++) {
    for (size_t j = 0; j < range.getSize(); j++) {
      sgpp::base::DataVector row(2);
      row.set(1, range.get(i));
      row.set(0, range.get(j));
      values.appendRow(row);
    }
  }
  // evaluate learned function at all points from values
  // and write result to csv file
  //output.open("SVM_ripley_errormonitor_impurity_fun_evals.csv");
  //output.open("SVM_ripley_periodicmonitor_impurity_fun_evals.csv");
  output.open("SVM_banana_errormonitor_combinedmeasure_fun_evals.csv");
  //output.open("SVM_banana_periodicmonitor_combinedmeasure_fun_evals.csv");
  size_t dim = values.getNcols();
  for (size_t i = 0; i < values.getNrows(); i++) {
    // Get next test sample x 
    sgpp::base::DataVector x(2);					
    values.getRow(i, x);
    //std::cout << x[0] << " , " << x[1] << std::endl;
    double res = svm->predictRaw(*grid, x, dim);
    output << res << ";" << std::endl;
  }  
  output.close();

}

double LearnerSVM::getAccuracy(sgpp::base::DataMatrix& testDataset,
                               const sgpp::base::DataVector& classesReference,
                               const double threshold) {
  // evaluate test dataset
  sgpp::base::DataVector classesComputed(testDataset.getNrows());
  predict(testDataset, classesComputed);

  return getAccuracy(classesComputed, classesReference, threshold);
}

double LearnerSVM::getAccuracy(const sgpp::base::DataVector& classesComputed,
                               const sgpp::base::DataVector& classesReference,
                               const double threshold) {
  double result = -1.0;

  if (classesComputed.getSize() != classesReference.getSize()) {
    throw base::application_exception(
        "LearnerSVM::getAccuracy: lengths of classes vectors do not match!");
  }

  size_t correct = 0;

  for (size_t i = 0; i < classesComputed.getSize(); i++) {
    if ((classesComputed.get(i) >= threshold && classesReference.get(i) >= 0.0) ||
       (classesComputed.get(i) < threshold && classesReference.get(i) < 0.0)) {
      correct++;
    }
  }

  result = static_cast<double>(correct) / static_cast<double>(classesComputed.getSize());

  return result;
}

void LearnerSVM::predict(sgpp::base::DataMatrix& testData,
                         sgpp::base::DataVector& computedLabels) {
  computedLabels.resize(testData.getNrows());
  size_t dim = testData.getNcols();  

  for (size_t i = 0; i < testData.getNrows(); i++) {
    // Get next test sample x 
    sgpp::base::DataVector x(dim);					
    testData.getRow(i, x);
    computedLabels.set(i, svm->predict(*grid, x, dim)); 
  }
}

double LearnerSVM::getError(sgpp::base::DataMatrix& data, sgpp::base::DataVector& labels, 
                            std::string errorType) {
  size_t numData = data.getNrows();
  size_t dim = data.getNcols();
  sgpp::base::DataVector error(numData);
  error.setAll(0.0);

  double res = -1.0;
  if (errorType == "MSE") {
    for (size_t i = 0; i < numData; i++) {
      sgpp::base::DataVector x(dim);					
      data.getRow(i, x);
      error.set(i, labels.get(i) - svm->predictRaw(*grid, x, dim));
    }
    // Loss (MSE)
    double sum = 0;
    for (size_t i = 0; i < numData; i++) {
      sum += error.get(i) * error.get(i);
    }
    res = (sum / (double)numData);
  }
  if (errorType == "Hinge") {
    for (size_t i = 0; i < numData; i++) {
      sgpp::base::DataVector x(dim);					
      data.getRow(i, x);
      error.set(i, std::max(0.0, 1.0 - labels.get(i) * svm->predictRaw(*grid, x, dim)));
    }
    // Loss (Hinge)
    double sum = 0;
    for (size_t i = 0; i < numData; i++) {
      sum += error.get(i);
    }
    res = (sum / (double)numData);
  }
  return res;
}

}  // namespace datadriven
}  // namespace sgpp


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
#include <sgpp/base/exception/data_exception.hpp>

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

LearnerSVM::LearnerSVM(sgpp::base::RegularGridConfiguration& gridConfig,
                       sgpp::base::AdpativityConfiguration& adaptConfig)                      
    : grid(nullptr),
      trainData(nullptr),
      trainLabels(nullptr),
      testData(nullptr),
      testLabels(nullptr),
      validData(nullptr),
      validLabels(nullptr),
      gridConfig(gridConfig),
      adaptivityConfig(adaptConfig),
      error(0.0) {}

LearnerSVM::~LearnerSVM() {}

void LearnerSVM::initialize(base::DataMatrix& pTrainData, base::DataVector& pTrainLabels,
                            base::DataMatrix& pTestData, base::DataVector& pTestLabels,
                            std::shared_ptr<sgpp::base::DataMatrix> pValidData,
                            std::shared_ptr<sgpp::base::DataVector> pValidLabels,
                            size_t budget) {
  trainData = std::make_shared<base::DataMatrix>(pTrainData);
  trainLabels = std::make_shared<base::DataVector>(pTrainLabels);
  testData = std::make_shared<base::DataMatrix>(pTestData);
  testLabels = std::make_shared<base::DataVector>(pTestLabels);
  validData = pValidData;
  validLabels = pValidLabels;
  gridConfig.dim_ = trainData->getNcols();
  grid = createRegularGrid();

  std::cout << "# initial grid size: " << grid->getSize() << std::endl;
  std::cout << "# SVS budget: " << budget << std::endl;
  
  // generate SVM
  svm = std::shared_ptr<PrimalDualSVM>(new PrimalDualSVM(
    grid->getSize(), trainData->getNcols(), budget, false));
}

std::shared_ptr<base::Grid> LearnerSVM::createRegularGrid() {
  // load grid
  std::unique_ptr<base::Grid> uGrid;
  if (gridConfig.type_ == base::GridType::ModLinear) {
    uGrid = base::Grid::createModLinearGrid(gridConfig.dim_);
  } 
  else {
    throw base::application_exception(
      "LearnerSVM::createRegularGrid : grid type is not supported");
  }

  uGrid->getGenerator().regular(gridConfig.level_);

  // move the grid to be shared
  std::shared_ptr<base::Grid> sGrid{std::move(uGrid)};

  return sGrid;
}

void LearnerSVM::train(size_t maxDataPasses, double lambda, double betaRef,
                       std::string refType, std::string refMonitor, 
                       size_t refPeriod, double errorDeclineThreshold,
                       size_t errorDeclineBufferSize, size_t minRefInterval) {

  GridStorage& gridStorage = grid->getStorage();

  size_t dim = trainData->getNcols();
  // initialize counter for dataset passes
  size_t cntDataPasses = 0;
  // iteration counter for computing learning parameter
  size_t t = 0;
 
  // initialize learning parameter 
  double eta = 0.0;
  // varible to store eta*classLabel
  double beta = 0.0;
  // raw prediction value
  double rawPrediction = 0.0;
  
  // refinement variables
  // counter for performed refinements
  size_t refCnt = 0;
  bool doRefine = false; // set true by monitor to trigger refinement
  double currentValidError = 0.0;
  double currentTrainError = 0.0;
  std::shared_ptr<ConvergenceMonitor> monitor(new ConvergenceMonitor(
    errorDeclineThreshold, errorDeclineBufferSize, minRefInterval));

  // auxiliary variable for accuracy (error) measurement
  double acc = getAccuracy(*testData, *testLabels, 0.0);
  avgErrors.append(1.0 - acc); 
 
  // counts total number of processed data points
  size_t processedPoints = 0;

  // learn the data
  while (cntDataPasses < maxDataPasses) {		
    for (size_t i = 0; i < trainData->getNrows(); i++) {
      // Get next training sample x and its label y
      sgpp::base::DataVector x(dim);					
      trainData->getRow(i, x);
      double y = trainLabels->get(i);                           
	  
      t += 1; 
      // compute next learning rate
      eta = 1.0 / (lambda * static_cast<double>(t));
      // update model with given data sample
      rawPrediction = svm->predictRaw(*grid, x, dim);
      svm->multiply(1.0 - lambda * eta);
      if (rawPrediction * y < 1.0) {
        beta = eta * y;
	svm->add(*grid, x, beta, dim);
      }	

      // check if refinement should be performed
      if (refMonitor == "convergence") {
        //convergence monitor
        if (refCnt < adaptivityConfig.numRefinements_) {
          if (validData == nullptr) {
            throw base::data_exception(
	      "No validation data for checking convergence provided!");
          }
          currentValidError = getError(*validData, *validLabels, "Hinge");
          currentTrainError = getError(*trainData, *trainLabels, "Hinge");
     
          monitor->pushToBuffer(currentValidError,currentTrainError);

          if (monitor->nextRefCnt > 0) {
            monitor->nextRefCnt--;
          }
          if (monitor->nextRefCnt == 0){
            doRefine = monitor->checkConvergence();
          }
        }
      }
      else if (refMonitor == "periodic") {
        // periodic monitor
        if ( (refCnt < adaptivityConfig.numRefinements_) 
        &&   (processedPoints > 0) 
        &&   ((processedPoints+1) % refPeriod == 0) ) {
          doRefine = true;
        }
      }

      processedPoints ++;

      if (doRefine) {
        //acc = getAccuracy(*testData, *testLabels, 0.0);
        //avgErrors.append(1.0 - acc);

        std::cout << "refinement at iteration: " << processedPoints << std::endl;
        HashRefinement refinement;

        if (refType == "combined-measure") {
          ForwardSelectorRefinement decorator(&refinement);
          // generate indicator
          ForwardSelectorRefinementIndicator indicator(
            *grid, *(svm->svs), *(svm->alphas), *(svm->w), 
            *(svm->w2), betaRef, adaptivityConfig.threshold_, 
            adaptivityConfig.noPoints_);
          // refine points according to indicator
          decorator.free_refine(gridStorage, indicator);
        }
        else if (refType == "impurity") {
          // compute current labels of support vectors 
          sgpp::base::DataVector svsClassesComputed((svm->svs)->getNrows());
          predict(*(svm->svs), svsClassesComputed);

          ImpurityRefinement decorator(&refinement);
          // generate indicator 
          ImpurityRefinementIndicator indicator(
            *grid, *(svm->svs), (svm->alphas).get(), (svm->w).get(),
            (svm->w2).get(), svsClassesComputed, adaptivityConfig.threshold_, 
            adaptivityConfig.noPoints_);
          // refine points according to indicator
          decorator.free_refine(gridStorage, indicator);
        }
        std::cout << "refinement step: " << refCnt+1 << std::endl;
        std::cout << "new grid size: " << grid->getSize() << std::endl;

        refCnt++;
        doRefine = false;
        if (refMonitor == "convergence") {
          monitor->nextRefCnt = monitor->minRefInterval;
        }
      }

      // save current error
      if (processedPoints % 10 == 0) {
        acc = getAccuracy(*testData, *testLabels, 0.0);
        avgErrors.append(1.0 - acc);
      }

    }
    cntDataPasses++;

  }
  std::cout << "final grid size: " << grid->getSize() << std::endl;
  double hinge = getError(*testData, *testLabels, "Hinge");
  std::cout << "hinge loss: " << hinge << std::endl;

  error = 1.0 - getAccuracy(*testData, *testLabels, 0.0);

}

void LearnerSVM::storeResults(sgpp::base::DataMatrix& testDataset) {

  sgpp::base::DataVector predictedLabels(testDataset.getNrows());
  predict(testDataset, predictedLabels);

  std::ofstream output;
  // write computed labels to csv file
  output.open("SVM_predicted_classes.csv");
  if (output.fail()) {
    std::cout << "failed to create csv file!" << std::endl;  
  }
  else {
    for (size_t i = 0; i < predictedLabels.getSize(); i++) {
      sgpp::base::DataVector x(2);					
      testDataset.getRow((size_t)i, x);
      output << x[0] << ";" << x[1] << ";" << predictedLabels[i] << std::endl;
    }
    output.close();
  }
  // write grid to csv file
  output.open("SVM_grid.csv");
  if (output.fail()) {
    std::cout << "failed to create csv file!" << std::endl;  
  }
  else {
    sgpp::base::GridStorage& storage = grid->getStorage();
    sgpp::base::GridStorage::grid_map_iterator end_iter = storage.end();
    for (sgpp::base::GridStorage::grid_map_iterator iter = storage.begin(); 
         iter != end_iter; iter++) { 
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

  // write function evaluations to csv file
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
  output.open("SVM_fun_evals.csv");
  size_t dim = values.getNcols();
  for (size_t i = 0; i < values.getNrows(); i++) {
    // get next test sample x 
    sgpp::base::DataVector x(2);					
    values.getRow(i, x);
    double res = svm->predictRaw(*grid, x, dim);
    output << res << ";" << std::endl;
  }  
  output.close();

}

double LearnerSVM::getAccuracy(sgpp::base::DataMatrix& testDataset,
                               const sgpp::base::DataVector& referenceLabels,
                               const double threshold) {
  // evaluate test dataset
  sgpp::base::DataVector predictedLabels(testDataset.getNrows());
  predict(testDataset, predictedLabels);

  return getAccuracy(referenceLabels, threshold, predictedLabels);
}

double LearnerSVM::getAccuracy(const sgpp::base::DataVector& referenceLabels,
                               const double threshold,
                               const sgpp::base::DataVector& predictedLabels) {
  double result = -1.0;

  if (predictedLabels.getSize() != referenceLabels.getSize()) {
    throw base::application_exception(
      "LearnerSVM::getAccuracy: lengths of classes vectors do not match!");
  }

  size_t correct = 0;

  for (size_t i = 0; i < predictedLabels.getSize(); i++) {
    if ((predictedLabels.get(i) >= threshold && referenceLabels.get(i) >= 0.0) ||
       (predictedLabels.get(i) < threshold && referenceLabels.get(i) < 0.0)) {
      correct++;
    }
  }

  result = static_cast<double>(correct) / static_cast<double>(predictedLabels.getSize());

  return result;
}

void LearnerSVM::predict(sgpp::base::DataMatrix& testData,
                         sgpp::base::DataVector& predictedLabels) {
  predictedLabels.resize(testData.getNrows());
  size_t dim = testData.getNcols();  

  for (size_t i = 0; i < testData.getNrows(); i++) {
    // get next test sample x 
    sgpp::base::DataVector x(dim);					
    testData.getRow(i, x);
    predictedLabels.set(i, svm->predict(*grid, x, dim)); 
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
    // loss (MSE)
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
    // loss (Hinge)
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


// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/LearnerSVM.hpp>

#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/ForwardSelectorRefinement.hpp>
#include <sgpp/base/grid/generation/functors/ForwardSelectorRefinementIndicator.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/ImpurityRefinement.hpp>
#include <sgpp/base/grid/generation/functors/ImpurityRefinementIndicator.hpp>

#include <sgpp/base/exception/application_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

////////////////////////////////////////////////////
#include "sgpp/datadriven/tools/ARFFTools.hpp"
////////////////////////////////////////////////////

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
      gridConfig(gridConfig),
      error(0.0) {}

LearnerSVM::~LearnerSVM() {}

// -----------------------------------------------------------------------

void LearnerSVM::initialize(base::DataMatrix& pTrainData, base::DataVector& pTrainLabels,
                            base::DataMatrix& pTestData, base::DataVector& pTestLabels) {
  trainData = std::make_shared<base::DataMatrix>(pTrainData);
  trainLabels = std::make_shared<base::DataVector>(pTrainLabels);
  testData = std::make_shared<base::DataMatrix>(pTestData);
  testLabels = std::make_shared<base::DataVector>(pTestLabels);
  gridConfig.dim_ = trainData->getNcols();
  grid = createRegularGrid();

  std::cout << "# grid points: " << grid->getSize() << std::endl;
  
  int B = 100; //ToDo: pass budget B (=max number of stored support vectors) as parameter
               // Ripley 100 / Banana 10000
  
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
  
  size_t maxIterations = 2; // Ripley 3 - Banana 1 // ToDo: pass as parameter
  size_t numIterations = 0;
  size_t t = 0;

  // refinement parameters - ToDo: pass to Learner as parameters (adaptivityConfig)
  size_t numRefSteps = 5;
  size_t refPeriod = 45; //Ripley 25/50 - Banana 100/
  size_t numPoints = 5;
  double threshold = 0.0;
  double betaRef = 2.0;
 
  double lambda = 0.4; //ToDo: pass as parameter / Ripley 0.0125 - Banana 0.00125
  double eta = 0.0;
  double beta = 0.0;
  double rawPrediction = 0.0;

  bool doRefine = false;      // set true by monitor to trigger refinement
  size_t refSteps = 0;

  // for error plotting
  sgpp::base::DataVector errors;
  //sgpp::base::DataVector grids;
  // for plotting only
  //double acc = getAccuracy(*trainData, *trainLabels, 0.0);
  double acc = getAccuracy(*testData, *testLabels, 0.0);
  errors.append(1.0 - acc);

  size_t steps = 0;
    while (numIterations < maxIterations) {			
      for (size_t i = 0; i < trainData->getNrows(); i++) {
        //std::cout << "i: " << i << std::endl;
        // Get next training sample x and its label y
        //int idx = getRandom(static_cast<int>(trainData->getNrows()) - 1);
        sgpp::base::DataVector x(dim);					
        //trainData->getRow((size_t)idx, x);
        trainData->getRow((size_t)i, x);
        double y = trainLabels->get(i);                           
	  
        t += 1; 
        eta = 1.0 / (lambda * static_cast<double>(t));

        rawPrediction = svm->predictRaw(*grid, x, dim);
        svm->multiply(1.0 - lambda * eta);
        if (rawPrediction * y < 1.0) {
          beta = eta * y;
	  svm->add(*grid, x, beta, dim);
        }	  
      
        //if ( (i > 0) && ((i+1) % refPeriod == 0) ) {
        if ( (steps > 0) && ((steps+1) % refPeriod == 0) ) {
          doRefine = true;
        }

        if ( (refSteps < numRefSteps) && doRefine ) {
          // compute current labels of support vectors 
          // (needed for impurity refinement)
          sgpp::base::DataVector svsClassesComputed((svm->svs)->getNrows());
          predict(*(svm->svs), svsClassesComputed);

          HashRefinement refinement;
          //ForwardSelectorRefinement decorator(&refinement);
          ImpurityRefinement decorator(&refinement);

          // refine points according to indicator
          //ForwardSelectorRefinementIndicator indicator(*grid, *(svm->svs), *(svm->alphas), *(svm->w),
          //                                             *(svm->w2), betaRef, threshold, numPoints);        
          ImpurityRefinementIndicator indicator(*grid, *(svm->svs), *(svm->alphas), *(svm->w),
                                                *(svm->w2), svsClassesComputed, threshold, numPoints);

          decorator.free_refine(gridStorage, indicator);

          std::cout << "refinement step: " << refSteps+1 << std::endl;
          //std::cout << "size of w :" << svm->w->getSize() << std::endl;
          std::cout << "new grid size: " << grid->getSize() << std::endl;

          refSteps++;
          doRefine = false;
        }
        // for plotting only
        //acc = getAccuracy(*trainData, *trainLabels, 0.0);
        acc = getAccuracy(*testData, *testLabels, 0.0);
        errors.append(1.0 - acc);

        steps++;
      }
      //refSteps = 0;
      numIterations++;

    }
    //write error evaluation to .csv
    std::ofstream output;
    output.open("SVM_ripley_err_rate_impurity_train_"+std::to_string(dataNum)+".csv");
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

//void LearnerSVM::printResults()

double LearnerSVM::getAccuracy(sgpp::base::DataMatrix& testDataset,
                               const sgpp::base::DataVector& classesReference,
                               const double threshold) {
  // evaluate test dataset
  sgpp::base::DataVector classesComputed(testDataset.getNrows());
  predict(testDataset, classesComputed);

  //write computed classes to .csv
  /*std::ofstream output;
  //output.open("banana_impurity_predicted_train_1.csv");
  //output.open("banana_measure_predicted_train_7.csv");
  output.open("ripley_impurity_predicted_train_.csv");
  //output.open("ripley_measure_predicted_train_1.csv");
  if (output.fail()) {
    std::cout << "failed to create .csv file!" << std::endl;  
  }
  else {
    for (size_t i = 0; i < classesComputed.getSize(); i++) {
      sgpp::base::DataVector x(2);					
      testDataset.getRow((size_t)i, x);
      output << x[0] << ";" << x[1] << ";" << classesComputed[i] << std::endl;
    }
    output.close();
  }
  //write grid to .csv
  //output.open("banana_impurity_grid_train_1.csv");
  //output.open("banana_measure_grid_train_7.csv");
  output.open("ripley_impurity_grid_train_.csv");
  //output.open("ripley_measure_grid_train_1.csv");
  if (output.fail()) {
    std::cout << "failed to create .csv file!" << std::endl;  
  }
  else {
    GridStorage& storage = grid->getStorage();
    GridStorage::grid_map_iterator end_iter = storage.end();
    for (GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter; iter++) { 
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
  }*/

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
    testData.getRow((size_t)i, x);
    computedLabels.set(i, svm->predict(*grid, x, dim)); 
  }
}

}  // namespace datadriven
}  // namespace sgpp


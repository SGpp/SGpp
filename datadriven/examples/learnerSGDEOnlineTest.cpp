// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <string>
#include <random>

#include "sgpp/globaldef.hpp"
#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/datadriven/application/LearnerSGDE.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/datadriven/application/RegularizationConfiguration.hpp"
//#include "sgpp/datadriven/application/GaussianKDE.hpp"
//#include "sgpp/datadriven/DatadrivenOpFactory.hpp"

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridStorage;

/*
int main(int argc, char** argv) {
  size_t totalSets = 10; //10
  double avgError = 0.0;
  for (size_t numSets = 0; numSets < totalSets; numSets++) {
    //std::string filename = "/media/sf_Downloads/MA/ripley/simple/fixed_train_seed42/ripley_train_"+std::to_string(numSets+1)+".arff";
    std::string filename = "/media/sf_Downloads/MA/banana/simple/fixed_train_seed42/banana_train_"+std::to_string(numSets+1)+".arff";

    // load training samples
    std::cout << "# loading file: " << filename << std::endl;
    sgpp::datadriven::Dataset trainDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
    sgpp::base::DataMatrix& trainData = trainDataset.getData();

    //normalize 
    //trainData.normalizeDimension(0);
    //trainData.normalizeDimension(1);

    // extract train classes
    sgpp::base::DataVector& trainLabels = trainDataset.getTargets();

    //filename = "/media/sf_Downloads/MA/ripley/simple/fixed_train_seed42/ripley_test.arff";
    filename = "/media/sf_Downloads/MA/banana/simple/fixed_train_seed42/banana_test.arff";

    // load test samples
    std::cout << "# loading file: " << filename << std::endl;
    sgpp::datadriven::Dataset testDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
    sgpp::base::DataMatrix& testData = testDataset.getData();

    //normalize 
    //testData.normalizeDimension(0);
    //testData.normalizeDimension(1);

    // extract test classes
    sgpp::base::DataVector& testLabels = testDataset.getTargets();

    // configure grid
    std::cout << "# create grid config" << std::endl;
    sgpp::base::RegularGridConfiguration gridConfig;
    gridConfig.dim_ = trainDataset.getDimension();
    gridConfig.level_ = 4;
    gridConfig.type_ = sgpp::base::GridType::Linear;
    //gridConfig.type_ = sgpp::base::GridType::ModLinear;

    // configure adaptive refinement
    std::cout << "# create adaptive refinement config" << std::endl;
    sgpp::base::AdpativityConfiguration adaptConfig;
    adaptConfig.numRefinements_ = 28;
    adaptConfig.noPoints_ = 3;

    // configure solver
    std::cout << "# create solver config" << std::endl;
    sgpp::solver::SLESolverConfiguration solverConfig;
    solverConfig.type_ = sgpp::solver::SLESolverType::CG;
    solverConfig.maxIterations_ = 20; // 1000
    solverConfig.eps_ = 1e-10; // 1e-14
    solverConfig.threshold_ = 1e-10; // 1e-14

    // configure regularization
    std::cout << "# create regularization config" << std::endl;
    sgpp::datadriven::RegularizationConfiguration regularizationConfig;
    //regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Identity;
    regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Laplace;

    // configure learner
    std::cout << "# create learner config" << std::endl;
    sgpp::datadriven::CrossvalidationForRegularizationConfiguration crossvalidationConfig;
    crossvalidationConfig.enable_ = false;
    crossvalidationConfig.lambda_ = 0.0015;

    crossvalidationConfig.kfold_ = 5; // 5
    crossvalidationConfig.lambdaStart_ = 1e-1;
    crossvalidationConfig.lambdaEnd_ = 1e-10;
    crossvalidationConfig.lambdaSteps_ = 5;
    crossvalidationConfig.logScale_ = true;
    crossvalidationConfig.shuffle_ = true;
    crossvalidationConfig.seed_ = 1234567;
    crossvalidationConfig.silent_ = true;

    std::cout << "# creating the learner" << std::endl;
    sgpp::datadriven::LearnerSGDE learner(gridConfig, adaptConfig, solverConfig, regularizationConfig,
                                        crossvalidationConfig);
    learner.initialize(trainData);
    
    bool usePrior = true;

    std::cout << "# start to train the learner" << std::endl;
    learner.trainOnline(trainLabels, testData, testLabels, nullptr, nullptr, numSets+1, usePrior);

    std::cout << "# finished training" << std::endl;

    // test learner
    double accTrain = learner.getAccuracy(trainData, trainLabels, 0.0);
    std::cout << "Acc (train): " << accTrain << std::endl;
    double accTest = learner.getAccuracy(testData, testLabels, 0.0);
    std::cout << "Acc (test): " << accTest << std::endl;

    //learner.storeResults(testData, testLabels, 0.0);

    avgError += learner.error;
  }
  avgError = avgError / static_cast<double>(totalSets);
  std::cout << "Average accuracy on test data: " << (1.0 - avgError) << std::endl;
}
*/

int main(int argc, char** argv) {
  size_t totalSets = 1; //10
  size_t totalFolds = 1;
  double avgError = 0.0;
  double avgErrorFolds = 0.0;
  for (size_t numSets = 0; numSets < totalSets; numSets++) {
    sgpp::base::DataVector avgErrorsFolds(22, 0.0);
    for (size_t numFolds = 0; numFolds < totalFolds; numFolds++) {
      //std::string filename = "/media/sf_Downloads/MA/ripley/5_fold/fixed_train_seed42/ripley_train_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      std::string filename = "/media/sf_Downloads/MA/banana/5_fold/fixed_train_seed42/banana_train_10_4.arff";
      //std::string filename = "/media/sf_Downloads/MA/banana/5_fold/fixed_train_seed42/banana_train_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      //std::string filename = "/media/sf_Downloads/MA/DR10/5_fold/fixed_train_seed42/50k/DR10_train_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      // load training samples
      std::cout << "# loading file: " << filename << std::endl;
      sgpp::datadriven::Dataset trainDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
      sgpp::base::DataMatrix& trainData = trainDataset.getData();
      // extract train classes
      sgpp::base::DataVector& trainLabels = trainDataset.getTargets();

      //filename = "/media/sf_Downloads/MA/ripley/5_fold/fixed_train_seed42/ripley_test.arff";
      filename = "/media/sf_Downloads/MA/banana/5_fold/fixed_train_seed42/banana_test.arff";
      //filename = "/media/sf_Downloads/MA/DR10/5_fold/fixed_train_seed42/50k/DR10_test.arff";
      // load test samples
      std::cout << "# loading file: " << filename << std::endl;
      sgpp::datadriven::Dataset testDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
      sgpp::base::DataMatrix& testData = testDataset.getData();
      // extract test classes
      sgpp::base::DataVector& testLabels = testDataset.getTargets();

      //filename = "/media/sf_Downloads/MA/ripley/5_fold/fixed_train_seed42/ripley_val_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      filename = "/media/sf_Downloads/MA/banana/5_fold/fixed_train_seed42/banana_val_10_4.arff";
      //filename = "/media/sf_Downloads/MA/banana/5_fold/fixed_train_seed42/banana_val_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      //filename = "/media/sf_Downloads/MA/DR10/5_fold/fixed_train_seed42/50k/DR10_val_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      // load validation samples
      std::cout << "# loading file: " << filename << std::endl;
      sgpp::datadriven::Dataset valDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
      sgpp::base::DataMatrix* validData = &valDataset.getData();
      // extract validation classes
      sgpp::base::DataVector* validLabels = &valDataset.getTargets();

      // configure grid
      std::cout << "# create grid config" << std::endl;
      sgpp::base::RegularGridConfiguration gridConfig;
      gridConfig.dim_ = trainDataset.getDimension();
      gridConfig.level_ = 3;
      gridConfig.type_ = sgpp::base::GridType::Linear;
      //gridConfig.type_ = sgpp::base::GridType::ModLinear;

      // configure adaptive refinement
      std::cout << "# create adaptive refinement config" << std::endl;
      sgpp::base::AdpativityConfiguration adaptConfig;
      adaptConfig.numRefinements_ = 20; 
      adaptConfig.noPoints_ = 10; 

      // configure solver
      std::cout << "# create solver config" << std::endl;
      sgpp::solver::SLESolverConfiguration solverConfig;
      solverConfig.type_ = sgpp::solver::SLESolverType::CG;
      solverConfig.maxIterations_ = 20; //1000
      solverConfig.eps_ = 1e-10; //1e-14
      solverConfig.threshold_ = 1e-10; //1e-14

      // configure regularization
      std::cout << "# create regularization config" << std::endl;
      sgpp::datadriven::RegularizationConfiguration regularizationConfig;
      //regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Identity;
      regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Laplace;

      // configure learner
      std::cout << "# create learner config" << std::endl;
      sgpp::datadriven::CrossvalidationForRegularizationConfiguration crossvalidationConfig;
      crossvalidationConfig.lambda_ = 0.0015;
      crossvalidationConfig.enable_ = false;

      crossvalidationConfig.kfold_ = 5; // 5
      crossvalidationConfig.lambdaStart_ = 1e-1;
      crossvalidationConfig.lambdaEnd_ = 1e-10;
      crossvalidationConfig.lambdaSteps_ = 5;
      crossvalidationConfig.logScale_ = true;
      crossvalidationConfig.shuffle_ = true;
      crossvalidationConfig.seed_ = 1234567;
      crossvalidationConfig.silent_ = true;

      std::cout << "# creating the learner" << std::endl;
      sgpp::datadriven::LearnerSGDE learner(gridConfig, adaptConfig, solverConfig, regularizationConfig,
                                        crossvalidationConfig);
      learner.initialize(trainData);
    
      bool usePrior = true;

      std::cout << "# start to train the learner" << std::endl;
      learner.trainOnline(trainLabels, testData, testLabels, validData, validLabels, numSets+1, usePrior);

      std::cout << "# finished training" << std::endl;

      // test learner
      double accTrain = learner.getAccuracy(trainData, trainLabels, 0.0);
      std::cout << "Acc (train): " << accTrain << std::endl;
      double accTest = learner.getAccuracy(testData, testLabels, 0.0);
      std::cout << "Acc (test): " << accTest << std::endl;

      //learner.storeResults(testData, testLabels, 0.0);

      avgErrorFolds += learner.error;
      avgErrorsFolds.add(learner.avgErrors);
    }
    avgErrorFolds = avgErrorFolds / static_cast<double>(totalFolds);
    std::cout << "Average accuracy on test data (set "+std::to_string(numSets+1)+"): " << (1.0 - avgErrorFolds) << std::endl;
    avgError += avgErrorFolds;
    avgErrorFolds = 0.0;

    avgErrorsFolds.mult(1.0/static_cast<double>(totalFolds));
    //write error evaluation to .csv
    std::ofstream output;
    //output.open("SGDE_ripley_5_fold_error_monitor_err_rate_zerocrossing_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SGDE_ripley_5_fold_error_monitor_err_rate_databased_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SGDE_ripley_5_fold_error_monitor_err_rate_norefinement_l3grid_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SGDE_banana_5_fold_error_monitor_err_rate_zerocrossing_train_l3startgrid_15epochs_"+std::to_string(numSets+1)+".csv");
    //output.open("SGDE_banana_5_fold_periodic_monitor_err_rate_zerocrossing_train_l3startgrid_15epochs_"+std::to_string(numSets+1)+".csv");
    //output.open("SGDE_banana_5_fold_err_rate_norefinement_train_l3startgrid_15epochs_"+std::to_string(numSets+1)+".csv");
    //output.open("SGDE_banana_5_fold_error_monitor_err_rate_databased_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SGDE_banana_5_fold_error_monitor_err_rate_norefinement_train_l3startgrid_10epochs_"+std::to_string(numSets+1)+".csv");
    //output.open("SGDE_banana_err_rate_zerocrossing_train_l3start"+std::to_string(numSets+1)+".csv");
    //output.open("SGDE_banana_err_rate_databased_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SGDE_DR10_5_fold_periodic_monitor_err_rate_zerocrossing_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SGDE_DR10_5_fold_error_monitor_err_rate_zerocrossing_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SGDE_DR10_5_fold_error_monitor_err_rate_zerocrossing_train_l5"+std::to_string(numSets+1)+".csv");
    //output.open("SGDE_DR10_5_fold_err_rate_norefinement_train_"+std::to_string(numSets+1)+".csv");
    output.open("SGDE_banana_5_fold_error_monitor_err_rate_zerocrossing_tocmpwithgridgrowth_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SGDE_banana_5_fold_error_monitor_err_rate_databased_tocmpwithgridgrowth_train_"+std::to_string(numSets+1)+".csv");
    if (output.fail()) {
      std::cout << "failed to create .csv file!" << std::endl;  
    }
    else {
      for (size_t i = 0; i < avgErrorsFolds.getSize(); i++) {					
        output << avgErrorsFolds.get(i) << ";" << std::endl;
      }
      output.close();
    }
  }

  avgError = avgError / static_cast<double>(totalSets);
  std::cout << "Average accuracy on test data: " << (1.0 - avgError) << std::endl;
}



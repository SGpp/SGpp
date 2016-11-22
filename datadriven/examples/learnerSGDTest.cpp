// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <string>

#include "sgpp/globaldef.hpp"
#include "sgpp/base/datatypes/DataVector.hpp"
#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/datadriven/application/LearnerSGD.hpp"
#include "sgpp/base/grid/Grid.hpp"

/*
int main(int argc, char** argv) {
  size_t totalSets = 10;
  double avgError = 0.0;
  for (size_t numSets = 0; numSets < totalSets; numSets++) {
    //std::string filename = "/media/sf_Downloads/MA/ripley/simple/fixed_train_seed42/ripley_train_"+std::to_string(numSets+1)+".arff";
    std::string filename = "/media/sf_Downloads/MA/banana/simple/fixed_train_seed42/banana_train_"+std::to_string(numSets+1)+".arff";

    // load training samples
    std::cout << "# loading file: " << filename << std::endl;
    sgpp::datadriven::Dataset trainDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
    sgpp::base::DataMatrix& trainData = trainDataset.getData();  

    // extract train classes
    sgpp::base::DataVector& trainLabels = trainDataset.getTargets();
  
    //filename = "/media/sf_Downloads/MA/ripley/simple/fixed_train_seed42/ripley_test.arff";
    filename = "/media/sf_Downloads/MA/banana/simple/fixed_train_seed42/banana_test.arff";

    // load test samples
    std::cout << "# loading file: " << filename << std::endl;
    sgpp::datadriven::Dataset testDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
    sgpp::base::DataMatrix& testData = testDataset.getData();

    // extract test classes
    sgpp::base::DataVector& testLabels = testDataset.getTargets();  

    // configure grid
    std::cout << "# creating grid config" << std::endl;
    sgpp::base::RegularGridConfiguration gridConfig;
    gridConfig.dim_ = trainDataset.getDimension();
    gridConfig.level_ = 4; // ripley 3 / banana 4
    //gridConfig.type_ = sgpp::base::GridType::Linear;
    gridConfig.type_ = sgpp::base::GridType::ModLinear;

    // configure adaptive refinement
    std::cout << "# create adaptive refinement config" << std::endl;
    sgpp::base::AdpativityConfiguration adaptConfig;
    adaptConfig.numRefinements_ = 25; 
    adaptConfig.noPoints_ = 5; 
    adaptConfig.threshold_ = 0.0;

    // specify additional parameters
    //size_t maxRuns = 1;
    double lambda = 0.00225; 
    double gamma = 0.11; 
    double errorDeclineThreshold = 1e-10; 
    size_t batchSize = 50;
    size_t bufferSize = 200; 
  
    // create SGD learner
    std::cout << "# creating the learner" << std::endl;
    sgpp::datadriven::LearnerSGD learner(gridConfig, adaptConfig);
  
    // initialize learner
    learner.initialize(trainData, trainLabels, testData, testLabels, nullptr, nullptr,
                       lambda, gamma, errorDeclineThreshold,
                       batchSize, bufferSize, false);

    // train learner
    std::cout << "# start to train the learner" << std::endl;
    learner.train(numSets+1);

    std::cout << "# finished training" << std::endl;

    //learner.storeResults(testData, testLabels, 0.0);

    // test learner
    double accTrain = learner.getAccuracy(trainData, trainLabels, 0.0);
    std::cout << "Acc (train): " << accTrain << std::endl;
    double accTest = learner.getAccuracy(testData, testLabels, 0.0);
    std::cout << "Acc (test): " << accTest << std::endl;

    avgError += learner.error;
  }

  avgError = avgError / static_cast<double>(totalSets);
  std::cout << "Average accuracy on test data: " << (1.0 - avgError) << std::endl;

}
*/

int main(int argc, char** argv) {
  size_t totalSets = 10; //10
  size_t totalFolds = 5; //5
  double avgError = 0.0;
  double avgErrorFolds = 0.0;
  for (size_t numSets = 0; numSets < totalSets; numSets++) {
    sgpp::base::DataVector avgErrorsFolds(721, 0.0); //Banana 721
    for (size_t numFolds = 0; numFolds < totalFolds; numFolds++) {
      //std::string filename = "/media/sf_Downloads/MA/ripley/5_fold/fixed_train_seed42/ripley_train_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      //std::string filename = "/media/sf_Downloads/MA/banana/5_fold/fixed_train_seed42/banana_train_2_3.arff";
      std::string filename = "/media/sf_Downloads/MA/banana/5_fold/fixed_train_seed42/banana_train_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      //std::string filename = "/media/sf_Downloads/MA/DR10/5_fold/fixed_train_seed42/50k/DR10_train_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
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
      //filename = "/media/sf_Downloads/MA/banana/5_fold/fixed_train_seed42/banana_val_2_3.arff";
      filename = "/media/sf_Downloads/MA/banana/5_fold/fixed_train_seed42/banana_val_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      //filename = "/media/sf_Downloads/MA/DR10/5_fold/fixed_train_seed42/50k/DR10_val_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      // load validation samples
      std::cout << "# loading file: " << filename << std::endl;
      sgpp::datadriven::Dataset valDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
      //sgpp::base::DataMatrix* valData = &valDataset.getData();
      std::shared_ptr<sgpp::base::DataMatrix> valData = std::make_shared<sgpp::base::DataMatrix>(valDataset.getData());
      // extract validation classes
      //sgpp::base::DataVector* valLabels = &valDataset.getTargets();
      std::shared_ptr<sgpp::base::DataVector> valLabels = std::make_shared<sgpp::base::DataVector>(valDataset.getTargets());

      // configure grid
      std::cout << "# creating grid config" << std::endl;
      sgpp::base::RegularGridConfiguration gridConfig;
      gridConfig.dim_ = trainDataset.getDimension();
      gridConfig.level_ = 3; 
      //gridConfig.type_ = sgpp::base::GridType::Linear;
      gridConfig.type_ = sgpp::base::GridType::ModLinear;

      // configure adaptive refinement
      std::cout << "# create adaptive refinement config" << std::endl;
      sgpp::base::AdpativityConfiguration adaptConfig;
      adaptConfig.numRefinements_ = 22; 
      adaptConfig.noPoints_ = 5;  
      adaptConfig.threshold_ = 0.0;

      // specify additional parameters
      //size_t maxRuns = 1;
      double lambda = 0.0015; 
      double gamma = 0.11; 
      double errorDeclineThreshold = 0.0005; 
      size_t batchSize = valData->getNrows(); 
      size_t bufferSize = 100; 
  
      // create SGD learner
      std::cout << "# creating the learner" << std::endl;
      sgpp::datadriven::LearnerSGD learner(gridConfig, adaptConfig);

      // initialize learner
      learner.initialize(trainData, trainLabels, testData, testLabels, valData, valLabels,
                         lambda, gamma, errorDeclineThreshold,
                         batchSize, bufferSize, true);
  
      // train learner
      std::cout << "# start to train the learner" << std::endl;
      learner.train(numSets+1);

      std::cout << "# finished training" << std::endl;

      //learner.storeResults(testData, testLabels, 0.0);

      // test learner
      double accTrain = learner.getAccuracy(trainData, trainLabels, 0.0);
      std::cout << "Acc (train): " << accTrain << std::endl;
      double accTest = learner.getAccuracy(testData, testLabels, 0.0);
      std::cout << "Acc (test): " << accTest << std::endl;

      avgErrorFolds += learner.error;
      //std::cout << avgErrorsFolds.getSize() << std::endl;
      //std::cout << learner.avgErrors.getSize() << std::endl;
      avgErrorsFolds.add(learner.avgErrors);
    }
    avgErrorFolds = avgErrorFolds / static_cast<double>(totalFolds);
    std::cout << "Average accuracy on test data (set "+std::to_string(numSets+1)+"): " << (1.0 - avgErrorFolds) << std::endl;
    avgError += avgErrorFolds;
    avgErrorFolds = 0.0;

    avgErrorsFolds.mult(1.0/static_cast<double>(totalFolds));
    //avgErrors.add(avgErrorsFolds);
    //write error evaluation to .csv
    std::ofstream output;
    //output.open("ASGD_ripley_5_fold_error_monitor_err_rate_predictive_train_"+std::to_string(numSets+1)+".csv");
    //output.open("ASGD_ripley_5_fold_error_monitor_err_rate_impurity_train_"+std::to_string(numSets+1)+".csv");
    //output.open("ASGD_ripley_5_fold_error_monitor_err_rate_norefinement_l3grid_train_"+std::to_string(numSets+1)+".csv");
    //output.open("ASGD_ripley_5_fold_error_monitor_Hinge_err_rate_predictive_train_"+std::to_string(numSets+1)+".csv");
    //output.open("ASGD_banana_5_fold_error_monitor_err_rate_predictive_train_l3startgrid_15epochs_"+std::to_string(numSets+1)+".csv");
    //output.open("ASGD_banana_5_fold_periodic_monitor_err_rate_predictive_train_l3startgrid_15epochs_"+std::to_string(numSets+1)+".csv");
    //output.open("ASGD_banana_5_fold_err_rate_norefinement_train_l3startgrid_15epochs_"+std::to_string(numSets+1)+".csv");
    //output.open("ASGD_banana_5_fold_error_monitor_err_rate_norefinement_l4grid_train_"+std::to_string(numSets+1)+".csv");
    //output.open("ASGD_banana_5_fold_error_monitor_err_rate_norefinement_l3grid_train_"+std::to_string(numSets+1)+".csv");
    //output.open("ASGD_banana_5_fold_error_monitor_Hinge_err_rate_predictive_train_"+std::to_string(numSets+1)+".csv");
    //output.open("ASGD_banana_err_rate_predictive_train_"+std::to_string(numSets+1)+".csv");
    //output.open("ASGD_banana_err_rate_impurity_train_"+std::to_string(numSets+1)+".csv");
    //output.open("ASGD_DR10_5_fold_error_monitor_err_rate_predictivve_train_"+std::to_string(numSets+1)+".csv");
    //output.open("ASGD_DR10_5_fold_error_monitor_err_rate_impurity_train_"+std::to_string(numSets+1)+".csv");
    //output.open("ASGD_DR10_5_fold_error_monitor_err_rate_impurity_train_l5"+std::to_string(numSets+1)+".csv");
    //output.open("ASGD_banana_5_fold_error_monitor_err_rate_predictive_tocmpwithgridgrowth_train_"+std::to_string(numSets+1)+".csv");
    //output.open("ASGD_banana_5_fold_error_monitor_err_rate_impurity_tocmpwithgridgrowth_train_"+std::to_string(numSets+1)+".csv");
    output.open("ASGD_banana_5_fold_error_monitor_err_rate_predictive_train_l3startgrid_15epochs_"+std::to_string(numSets+1)+".csv");
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






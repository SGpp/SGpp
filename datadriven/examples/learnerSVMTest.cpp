// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <string>

#include "sgpp/globaldef.hpp"
#include "sgpp/base/datatypes/DataVector.hpp"
#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/datadriven/application/LearnerSVM.hpp"
#include "sgpp/base/grid/Grid.hpp"

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

    // specify parameters
    //size_t maxRuns = 3;
    //ToDo:
  
    // create SVM learner
    std::cout << "# creating the learner" << std::endl;
    sgpp::datadriven::LearnerSVM learner(gridConfig);

    // initialize learner
    learner.initialize(trainData, trainLabels, testData, testLabels, nullptr, nullptr);

    // train learner
    std::cout << "# start to train the learner" << std::endl;
    learner.train(numSets+1);

    std::cout << "# finished training" << std::endl;

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
  size_t totalFolds = 5;
  double avgError = 0.0;
  double avgErrorFolds = 0.0;
  //size_t set = 0;
  for (size_t numSets = 0; numSets < totalSets; numSets++) {
    /*if (numSets == 0) {
      set = 0;
    }
    if (numSets == 1) {
      set = 2;
    }
    if (numSets == 2) {
      set = 4;
    }*/
    sgpp::base::DataVector avgErrorsFolds(721, 0.0); //banana 721 / DR10 12-17
    for (size_t numFolds = 0; numFolds < totalFolds; numFolds++) {
      //std::string filename = "/media/sf_Downloads/MA/ripley/5_fold/fixed_train_seed42/ripley_train_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      std::string filename = "/media/sf_Downloads/MA/banana/5_fold/fixed_train_seed42/banana_train_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      //std::string filename = "/media/sf_Downloads/MA/banana/5_fold/fixed_train_seed42/banana_train_10_1.arff";
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
      filename = "/media/sf_Downloads/MA/banana/5_fold/fixed_train_seed42/banana_val_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      //filename = "/media/sf_Downloads/MA/banana/5_fold/fixed_train_seed42/banana_val_10_1.arff";
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
      gridConfig.level_ = 3; //ripley,banana 3 / DR10 5
      //gridConfig.type_ = sgpp::base::GridType::Linear;
      gridConfig.type_ = sgpp::base::GridType::ModLinear;

      // specify parameters
      //size_t maxRuns = 3;
      //ToDo:
  
      // create SVM learner
      std::cout << "# creating the learner" << std::endl;
      sgpp::datadriven::LearnerSVM learner(gridConfig);

      // initialize learner
      learner.initialize(trainData, trainLabels, testData, testLabels, valData, valLabels);

      // train learner
      std::cout << "# start to train the learner" << std::endl;
      learner.train(numSets+1);

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
    //output.open("SVM_ripley_5_fold_error_monitor_err_rate_combinedmeasure_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SVM_ripley_5_fold_error_monitor_err_rate_impurity_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SVM_ripley_5_fold_error_monitor_err_rate_combinedmeasure_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SVM_ripley_5_fold_error_monitor_err_rate_norefinement_l3grid_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SVM_ripley_5_fold_error_monitor_Hinge_err_rate_impurity_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SVM_banana_5_fold_error_monitor_err_rate_combinedmeasure_train_l3startgrid_15epochs_"+std::to_string(numSets+1)+".csv");
    //output.open("SVM_banana_5_fold_error_monitor_err_rate_impurity_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SVM_banana_5_fold_error_monitor_err_rate_norefinement_l3grid_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SVM_banana_5_fold_error_monitor_Hinge_err_rate_combinedmeasure_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SVM_banana_5_fold_periodic_monitor_err_rate_combinedmeasure_train_l3startgrid_15epochs_"+std::to_string(numSets+1)+".csv");
    //output.open("SVM_banana_5_fold_err_rate_norefinement_train_l3startgrid_15epochs_"+std::to_string(numSets+1)+".csv");
    //output.open("SVM_banana_err_rate_impurity_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SVM_DR10_5_fold_periodic_monitor_err_rate_combinedmeasure_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SVM_DR10_5_fold_error_monitor_err_rate_combinedmeasure_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SVM_DR10_5_fold_error_monitor_err_rate_impurity_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SVM_DR10_5_fold_error_monitor_err_rate_impurity_train_l5"+std::to_string(numSets+1)+".csv");
    //output.open("SVM_DR10_5_fold_err_rate_norefinement_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SVM_banana_5_fold_error_monitor_err_rate_combinedmeasure_tocmpwithgridgrowth_train_"+std::to_string(numSets+1)+".csv");
    //output.open("SVM_banana_5_fold_error_monitor_err_rate_impurity_tocmpwithgridgrowth_train_"+std::to_string(numSets+1)+".csv");
    output.open("SVM_banana_5_fold_error_monitor_err_rate_combinedmeasure_train_l3startgrid_15epochs_"+std::to_string(numSets+1)+".csv");
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



// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/application/LearnerSGD.hpp>

#include <string>

int main() {
  size_t totalSets = 1; // 1-10 possible (10 differently ordered data sets)
  size_t totalFolds = 1; // set to 5 to perform 5-fold cv
  //double avgError = 0.0;
  //double avgErrorFolds = 0.0;
  for (size_t numSets = 0; numSets < totalSets; numSets++) {
    //sgpp::base::DataVector avgErrorsFolds(21, 0.0); // to compute average classification error
    for (size_t numFolds = 0; numFolds < totalFolds; numFolds++) {
      std::string filename = "../tests/data/ripley/5_fold/ripley_train_"
        +std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      //std::string filename = "../tests/data/banana/5_fold/banana_train_"
      //  +std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      //std::string filename = "../tests/data/SDSS_DR10/5_fold/DR10_train_"
      //  +std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      // load training samples
      std::cout << "# loading file: " << filename << std::endl;
      sgpp::datadriven::Dataset trainDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
      sgpp::base::DataMatrix& trainData = trainDataset.getData();  
      // extract train classes
      sgpp::base::DataVector& trainLabels = trainDataset.getTargets();
  
      filename = "../tests/data/ripley/5_fold/ripley_test.arff";
      //filename = "../tests/data/banana/5_fold/banana_test.arff";
      //filename = "../tests/data/SDSS_DR10/5_fold/DR10_test.arff";
      // load test samples
      std::cout << "# loading file: " << filename << std::endl;
      sgpp::datadriven::Dataset testDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
      sgpp::base::DataMatrix& testData = testDataset.getData();
      // extract test classes
      sgpp::base::DataVector& testLabels = testDataset.getTargets();  

      std::shared_ptr<sgpp::base::DataMatrix> validData = nullptr;
      std::shared_ptr<sgpp::base::DataVector> validLabels = nullptr;
      //if fixed validation data should be used (required for convergence monitor):
      filename = "../tests/data/ripley/5_fold/ripley_val_"
        +std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      //filename = "../tests/data/banana/5_fold/banana_val_"
      //  +std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      //filename = "../tests/data/SDSS_DR10/5_fold/DR10_val_"
      //  +std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      // load validation samples
      std::cout << "# loading file: " << filename << std::endl;
      sgpp::datadriven::Dataset valDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
      validData = std::make_shared<sgpp::base::DataMatrix>(valDataset.getData());
      // extract validation classes
      validLabels = std::make_shared<sgpp::base::DataVector>(valDataset.getTargets());

      // configure grid
      std::cout << "# creating grid config" << std::endl;
      sgpp::base::RegularGridConfiguration gridConfig;
      gridConfig.dim_ = trainDataset.getDimension();
      gridConfig.level_ = 3; 
      //gridConfig.type_ = sgpp::base::GridType::Linear;
      gridConfig.type_ = sgpp::base::GridType::ModLinear;

      // configure adaptive refinement
      std::cout << "# create adaptive refinement config" << std::endl;
      // possible refinement  monitors:
      // periodic monitor, convergence monitor
      std::string refMonitor;
      // select periodic monitor - perform refinements in fixed intervals
      // refMonitor = "periodic";
      size_t refPeriod = 50; // the refinement interval
      // select convergence monitor - perform refinements if algorithm has converged
      // (convergence measured with respect to MSE or Hinge loss observations)
      refMonitor = "convergence";
      double errorDeclineThreshold = 0.001; // the convergence threshold
      size_t errorDeclineBufferSize = 50; // number of error measurements which 
                                          // are considered for convergence check
      size_t minRefInterval = 0; // minimum number of iterations before next refinement 
                                 // is allowed to be performed
      std::cout << "Refinement monitor: " << refMonitor << std::endl;
      // possible refinement indicators:
      // predictive refinement, impurity-based refinement
      std::string refType;
      // select predictive refinement
      //refType = "predictive";
      // select impurity-based refinement
      refType = "impurity";
      std::cout << "Refinement type: " << refType << std::endl;
      sgpp::base::AdpativityConfiguration adaptConfig;
      adaptConfig.numRefinements_ = 2; 
      adaptConfig.noPoints_ = 5;  
      adaptConfig.threshold_ = 0.0;

      // additional parameters
      // specify max number of passes over traininig data set
      size_t maxDataPasses = 3;
      // regularization parameter
      double lambda = 0.3; 
      // initial learning rate
      double gamma = 0.5;  
      
      // specify number of data points to compute error contributions
      // for predictive refinement
      size_t batchSize = 50;
      if (validData != nullptr) {
        batchSize = validData->getNrows();
      }  
  
      // create SGD learner
      std::cout << "# creating the learner" << std::endl;
      sgpp::datadriven::LearnerSGD learner(gridConfig, adaptConfig);

      // initialize learner
      learner.initialize(trainData, trainLabels, testData, testLabels, validData, validLabels,
                         lambda, gamma, batchSize, true);
  
      // train learner
      std::cout << "# start to train the learner" << std::endl;
      learner.train(maxDataPasses, refType, refMonitor, 
                    refPeriod, errorDeclineThreshold,
                    errorDeclineBufferSize, minRefInterval);

      // store results (classified data, grid, function evaluations)
      //learner.storeResults(testData);

      // test learner
      double accTrain = learner.getAccuracy(trainData, trainLabels, 0.0);
      std::cout << "Acc (train): " << accTrain << std::endl;
      double accTest = learner.getAccuracy(testData, testLabels, 0.0);
      std::cout << "Acc (test): " << accTest << std::endl;

      //avgErrorFolds += learner.error;
      //avgErrorsFolds.add(learner.avgErrors);
    }
    /*avgErrorFolds = avgErrorFolds / static_cast<double>(totalFolds);
    std::cout << "Average accuracy on test data (set "+std::to_string(numSets+1)+"): " << (1.0 - avgErrorFolds) << std::endl;
    avgError += avgErrorFolds;
    avgErrorFolds = 0.0;

    avgErrorsFolds.mult(1.0/static_cast<double>(totalFolds));
    //avgErrors.add(avgErrorsFolds);*/

    //write error evaluation to csv file
    /*std::ofstream output;
    output.open("ASGD_avg_classification_error_"+std::to_string(numSets+1)+".csv");
    if (output.fail()) {
      std::cout << "failed to create .csv file!" << std::endl;  
    }
    else {
      for (size_t i = 0; i < avgErrorsFolds.getSize(); i++) {					
        output << avgErrorsFolds.get(i) << ";" << std::endl;
      }
      output.close();
    }*/
  }
  //avgError = avgError / static_cast<double>(totalSets);
  //std::cout << "Average accuracy on test data: " << (1.0 - avgError) << std::endl;

}






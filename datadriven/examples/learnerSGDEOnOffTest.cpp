// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/application/LearnerSGDEOnOff.hpp>
#include <sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp>

#include <string>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

int main() {
  size_t totalSets = 1; 
  size_t totalFolds = 1; 
  //double avgError = 0.0;
  //double avgErrorFolds = 0.0;
  for (size_t numSets = 0; numSets < totalSets; numSets++) {
    //sgpp::base::DataVector avgErrorsFolds(21, 0.0); //to compute average classification error
    for (size_t numFolds = 0; numFolds < totalFolds; numFolds++) {
      std::string filename = "/media/sf_Downloads/MA/ripley/5_fold/fixed_train_seed42/ripley_train_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      //std::string filename = "/media/sf_Downloads/MA/banana/5_fold/fixed_train_seed42/banana_train_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      //std::string filename = "/media/sf_Downloads/MA/DR10/5_fold/fixed_train_seed42/50k/DR10_train_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      // load training samples
      std::cout << "# loading file: " << filename << std::endl;
      sgpp::datadriven::Dataset trainDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
      sgpp::base::DataMatrix& trainData = trainDataset.getData();
      // extract train classes
      sgpp::base::DataVector& trainLabels = trainDataset.getTargets();

      filename = "/media/sf_Downloads/MA/ripley/5_fold/fixed_train_seed42/ripley_test.arff";
      //filename = "/media/sf_Downloads/MA/banana/5_fold/fixed_train_seed42/banana_test.arff";
      //filename = "/media/sf_Downloads/MA/DR10/5_fold/fixed_train_seed42/50k/DR10_test.arff";
      // load test samples
      std::cout << "# loading file: " << filename << std::endl;
      sgpp::datadriven::Dataset testDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
      sgpp::base::DataMatrix& testData = testDataset.getData();
      // extract test classes
      sgpp::base::DataVector& testLabels = testDataset.getTargets();

      sgpp::base::DataMatrix* validData = nullptr;
      sgpp::base::DataVector* validLabels = nullptr;
      //if fixed validation data should be used (required for convergence monitor):
      filename = "/media/sf_Downloads/MA/ripley/5_fold/fixed_train_seed42/ripley_val_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      //filename = "/media/sf_Downloads/MA/banana/5_fold/fixed_train_seed42/banana_val_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      //filename = "/media/sf_Downloads/MA/DR10/5_fold/fixed_train_seed42/50k/DR10_val_"+std::to_string(numSets+1)+"_"+std::to_string(numFolds+1)+".arff";
      // load validation samples
      std::cout << "# loading file: " << filename << std::endl;
      sgpp::datadriven::Dataset valDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
      validData = &valDataset.getData();
      // extract validation classes
      validLabels = &valDataset.getTargets();

      // set number of classes
      unsigned int classNum = 2;
 
      // set class labels
      sgpp::base::DataVector classLabels(classNum);
      classLabels[0] = -1;
      classLabels[1] = 1;  

      // configure grid
      std::cout << "# create grid config" << std::endl;
      sgpp::base::RegularGridConfiguration gridConfig;
      gridConfig.dim_ = trainDataset.getDimension();
      gridConfig.level_ = 3;
      gridConfig.type_ = sgpp::base::GridType::Linear;
      //gridConfig.type_ = sgpp::base::GridType::ModLinear;

      // configure regularization
      std::cout << "# create regularization config" << std::endl;
      sgpp::datadriven::RegularizationConfiguration regularizationConfig;
      regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Identity;

      //Define decomposition type
      DBMatDecompostionType dt;
      std::string decompType;
      //choose "LU decomposition"
      //dt = DBMatDecompLU;
      //decompType = "LU decomposition";
      //choose"Eigen decomposition"
      //dt = DBMatDecompEigen;
      //decompType = "Eigen decomposition";
      //choose "Cholesky decomposition"
      dt = DBMatDecompChol;
      decompType = "Cholesky decomposition";
      std::cout << "Decomposition type: " << decompType << std::endl;

      //if cholesky choosen -> configure adaptive refinement
      std::cout << "# create adaptive refinement configuration" << std::endl;
      //possible refinement  monitors:
      //periodic monitor, convergence monitor
      std::string refMonitor;
      //select periodic monitor - perform refinements in fixed intervals
      //refMonitor = "periodic";
      size_t refPeriod = 50; //the refinement interval
      //select convergence monitor - perform refinements if algorithm has converged
      //(convergence measured with respect to changes of the classification accuracy)
      refMonitor = "convergence";
      double accDeclineThreshold = 0.0005; //the convergence threshold
      size_t accDeclineBufferSize = 50; //number of accuracy measurements which 
                                         //are considered for convergence check
      size_t minRefInterval = 0; //minimum number of iterations before next refinement 
                                 //is allowed to be performed
      std::cout << "Refinement monitor: " << refMonitor << std::endl;
      //possible refinement indicators:
      //surplus refinement, data-based refinement, zero-crossings-based refinement
      std::string refType;
      //select surplus refinement
      //refType = "surplus";
      //select data-based refinement
      //refType = "data";
      //select zero-crossings-based refinement
      refType = "zero";
      std::cout << "Refinement type: " << refType << std::endl;
      sgpp::base::AdpativityConfiguration adaptConfig;
      adaptConfig.numRefinements_ = 2;
      adaptConfig.noPoints_ = 5;
      adaptConfig.threshold_ = 0.0; //only required for surplus refinement

      //initial lambda
      double lambda = 0.1; 

      double beta = 0.0; 

      sgpp::datadriven::DBMatDensityConfiguration dconf(&gridConfig, &adaptConfig, 
                                                        regularizationConfig.regType_, 
                                                        lambda, dt);

      bool usePrior = false; 

      //create learner
      std::cout << "# create learner" << std::endl;
      sgpp::datadriven::LearnerSGDEOnOff learner(dconf, trainData, trainLabels, 
                                                 testData, testLabels, validData, validLabels,
                                                 classLabels, classNum, lambda, usePrior, 
                                                 beta);

      //enable cv during learning
      bool enableCv = false;            
      //set cv configuration if cv enabled
      unsigned int nextCvStep = 50;
      double cvLambdaStart = 1e-1; 
      double cvLambdaEnd = 1e-10; 
      int cvLambdaSteps = 10;
      bool cvLogScale = true;
      sgpp::base::DataMatrix* cvTestData = &testData;
      sgpp::base::DataMatrix* cvTestDataRes = nullptr; //needed?
      learner.setCrossValidationParameters(cvLambdaSteps, cvLambdaStart, cvLambdaEnd, 
                                           cvTestData, cvTestDataRes, cvLogScale);

      //specify batch size
      size_t batchSize = 1;  
      //specify max number of passes over train data set
      size_t maxDataPasses = 3;

      std::cout << "# start to train the learner" << std::endl;
      learner.train(batchSize, maxDataPasses, refType, 
                    refMonitor, refPeriod, accDeclineThreshold,
                    accDeclineBufferSize, minRefInterval,
                    enableCv, nextCvStep);

      //compute accuracy
      double acc = learner.getAccuracy();
      std::cout << "# accuracy (test data): " << acc << std::endl;

      //store results (classified data, grids, density functions)
      //learner.storeResults();

      //avgErrorFolds += learner.error;
      //avgErrorsFolds.add(learner.avgErrors);
    }
    /*avgErrorFolds = avgErrorFolds / static_cast<double>(totalFolds);
    std::cout << "Average accuracy on test data (set "+std::to_string(numSets+1)+"): " << (1.0 - avgErrorFolds) << std::endl;
    avgError += avgErrorFolds;
    avgErrorFolds = 0.0;
    avgErrorsFolds.mult(1.0/static_cast<double>(totalFolds));*/

    //write error evaluation to csv-file
    /*std::ofstream output;
    output.open("SGDEOnOff_avg_classification_error_"+std::to_string(numSets+1)+".csv");
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
 
}


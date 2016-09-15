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


int main(int argc, char** argv) {
  size_t numSets = 10;
  double avgError = 0.0;
  for (size_t numSets = 0; numSets < 10; numSets++) {
    std::string filename = "../tests/data/ripley/ripley_train_"+std::to_string(numSets+1)+".arff";
    //std::string filename = "../tests/data/banana/banana_train_"+std::to_string(numSets+1)+".arff";

    //std::string filename = "../tests/data/ripley/ripley_train_1.arff";
    //std::string filename = "../tests/data/ripleyGarcke.train.arff";
    //std::string filename = "../tests/data/banana/banana_train_1.arff";
    //std::string filename = "../tests/data/banana.arff";
    // load training samples
    std::cout << "# loading file: " << filename << std::endl;
    sgpp::datadriven::Dataset trainDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
    sgpp::base::DataMatrix& trainData = trainDataset.getData();  

    //normalize features
    trainData.normalizeDimension(0);
    trainData.normalizeDimension(1);

    // extract train classes
    sgpp::base::DataVector& trainLabels = trainDataset.getTargets();
  
    filename = "../tests/data/ripley/ripley_test.arff";
    //filename = "../tests/data/ripleyGarcke.test.arff";
    //filename = "../tests/data/banana/banana_train_0.arff";
    //filename = "../tests/data/banana.arff";
    // load test samples
    std::cout << "# loading file: " << filename << std::endl;
    sgpp::datadriven::Dataset testDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
    sgpp::base::DataMatrix& testData = testDataset.getData();

    //normalize features
    testData.normalizeDimension(0);
    testData.normalizeDimension(1);

    // extract test classes
    sgpp::base::DataVector& testLabels = testDataset.getTargets();  

    // configure grid
    std::cout << "# creating grid config" << std::endl;
    sgpp::base::RegularGridConfiguration gridConfig;
    gridConfig.dim_ = trainDataset.getDimension();
    gridConfig.level_ = 2; 
    //gridConfig.type_ = sgpp::base::GridType::Linear;
    gridConfig.type_ = sgpp::base::GridType::ModLinear;

    // configure adaptive refinement
    std::cout << "# create adaptive refinement config" << std::endl;
    sgpp::base::AdpativityConfiguration adaptConfig;
    adaptConfig.numRefinements_ = 8; //banana 20 | ripley 7
    adaptConfig.noPoints_ = 7; //banana 5 | ripley 5
    adaptConfig.threshold_ = 0.0;

    // specify additional parameters
    //size_t maxRuns = 1;
    double lambda = 0.3; //1e-4 //banana 0.001 | ripley 0.005
    double gamma = 0.005; // banana 0.05 | ripley 0.005
    double errorDeclineThreshold = 0.000001; //1e-10
    size_t batchSize = 50;
    size_t bufferSize = 75; //banana 200 | ripley 75
  
    // create SVM learner
    std::cout << "# creating the learner" << std::endl;
    sgpp::datadriven::LearnerSGD learner(gridConfig, adaptConfig);

    // initialize learner
    learner.initialize(trainData, trainLabels, testData, testLabels,
                       lambda, gamma, errorDeclineThreshold,
                       batchSize, bufferSize);

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

  avgError = avgError / static_cast<double>(numSets);
  std::cout << "Average accuracy on test data: " << (1.0 - avgError) << std::endl;

}





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
//#include "sgpp/datadriven/DatadrivenOpFactory.hpp"


int main(int argc, char** argv) {
  std::string filename = "../tests/data/ripleyGarcke.train.arff";
  //std::string filename = "../tests/data/banana.train.arff";
  //std::string filename = "../tests/data/banana.arff";
  // load training samples
  std::cout << "# loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset trainDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
  sgpp::base::DataMatrix& trainData = trainDataset.getData();  

  //normalize - banana
  //trainData.normalizeDimension(0);
  //trainData.normalizeDimension(1);

  // extract train classes
  sgpp::base::DataVector& trainLabels = trainDataset.getTargets();

  filename = "../tests/data/ripleyGarcke.test.arff";
  //filename = "../tests/data/banana.test.arff";
  //filename = "../tests/data/banana.arff";
  // load test samples
  std::cout << "# loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset testDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
  sgpp::base::DataMatrix& testData = testDataset.getData();

  //normalize - banana
  //testData.normalizeDimension(0);
  //testData.normalizeDimension(1);

  // extract test classes
  sgpp::base::DataVector& testLabels = testDataset.getTargets();  

  // configure grid
  std::cout << "# creating grid config" << std::endl;
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = trainDataset.getDimension();
  gridConfig.level_ = 2; //ripley - fixed grid: l=4
  //gridConfig.type_ = sgpp::base::GridType::Linear;
  gridConfig.type_ = sgpp::base::GridType::ModLinear;

  // specify parameters
  //size_t maxIterations = 3;
  //ToDo:
  
  // create SVM learner
  std::cout << "# creating the learner" << std::endl;
  sgpp::datadriven::LearnerSVM learner(gridConfig);

  // initialize learner
  learner.initialize(trainData, trainLabels, testData, testLabels);

  // train learner
  std::cout << "# start training the learner" << std::endl;
  learner.train();

  std::cout << "# finished training" << std::endl;

  // test learner
  double accTrain = learner.getAccuracy(trainData, trainLabels, 0.0);
  std::cout << "Acc (train): " << accTrain << std::endl;
  double accTest = learner.getAccuracy(testData, testLabels, 0.0);
  std::cout << "Acc (test): " << accTest << std::endl;

}



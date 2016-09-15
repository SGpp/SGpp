// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/application/LearnerSGDEOnOff.hpp>
#include <sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp>

#include <ctime>
#include <string>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

int main(int argc, char** argv) {
  size_t totalSets = 10; // 10
  double avgError = 0.0;
  for (size_t numSets = 0; numSets < totalSets; numSets++) {
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

    // normalize data
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

    // normalize data
    testData.normalizeDimension(0);
    testData.normalizeDimension(1);

    // extract test classes
    sgpp::base::DataVector& testLabels = testDataset.getTargets();

    // set number of classes
    int classNum = 2;
 
    // set class labels
    /*double labels[classNum]; //correct??
    labels[0] = -1.0;
    labels[1] = 1.0;*/
    //ToDo: e.g. labels = [-1.0,1.0]
  

    // configure grid
    std::cout << "# create grid config" << std::endl;
    sgpp::base::RegularGridConfiguration gridConfig;
    gridConfig.dim_ = trainDataset.getDimension();
    gridConfig.level_ = 3;
    gridConfig.type_ = sgpp::base::GridType::Linear;

    // configure regularization
    std::cout << "# create regularization config" << std::endl;
    sgpp::datadriven::RegularizationConfiguration regularizationConfig;
    regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Identity;

    //Define decomposition type
    DBMatDecompostionType dt;
    std::string decompType;
    // "LU decomposition"
    //dt = DBMatDecompLU;
    //decompType = "LU decomposition";
    // "Eigen decomposition"
    //dt = DBMatDecompEigen;
    //decompType = "Eigen decomposition";
    //"Cholesky decomposition"
    dt = DBMatDecompChol;
    decompType = "Cholesky decomposition";
    std::cout << "Decomposition type: " << decompType << std::endl;

    // if cholesky choosen -> configure adaptive refinement
    std::cout << "# create adaptive refinement config" << std::endl;
    sgpp::base::AdpativityConfiguration adaptConfig;
    adaptConfig.numRefinements_ = 4;
    adaptConfig.noPoints_ = 6;
    adaptConfig.threshold_ = 0.0; //only required for surplus refinement

    // initial lambda
    double lambda = 2.5e-1; //1e-1 //ToDo: set to reasonable value

    double beta = 0.0; //ToDo: set to reasonable value

    sgpp::datadriven::DBMatDensityConfiguration dconf(&gridConfig, &adaptConfig, 
                                  regularizationConfig.regType_, lambda, dt);

    // cv configuration
    double cv_lambdaStart = 5e-1; //1e-1
    double cv_lambdaEnd = 2e-1; //1e-10
    int cv_lambdaSteps = 5;
    bool cv_logScale = true;

    bool usePrior = true; // false

    // create learner
    std::cout << "# create learner" << std::endl;
    LearnerSGDEOnOff learner(dconf, trainData, trainLabels, 
                             testData, testLabels,
                             classNum, lambda, usePrior, 
                             beta);

    //std::shared_ptr<sgpp::base::DataMatrix> cv_testData = std::make_shared<sgpp::base::DataMatrix>(testData);
    sgpp::base::DataMatrix* cv_testData = &testData;
    sgpp::base::DataMatrix* cv_testDataRes = nullptr; //ToDo: needed?
    learner.setCrossValidationParameters(cv_lambdaSteps, cv_lambdaStart, cv_lambdaEnd, 
                                         cv_testData, cv_testDataRes, cv_logScale);

    size_t batch_size = 1; // 1
    unsigned int next_cv_step = 100; //50
    std::cout << "# start to train the learner" << std::endl;
    learner.train(batch_size, next_cv_step, numSets+1);

    //compute accuracy
    double acc = learner.getAccuracy();
    std::cout << "# accuracy (test data): " << acc << std::endl;

    avgError += learner.error;
  }

  if (avgError > 0.0) {
    avgError = avgError / static_cast<double>(totalSets);
  }
  std::cout << "Average accuracy on test data: " << (1.0 - avgError) << std::endl;
}



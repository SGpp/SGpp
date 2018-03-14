// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/application/LearnerSGDEOnOff.hpp>
#include <sgpp/datadriven/configuration/DensityEstimationConfiguration.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>

#include <string>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

/**
 * This example shows how to perform offline/online-classification using sparse
 * grid density estimation on a interaction based grid
 * and matrix decomposition methods. It creates an
 * instance of LearnerSGDEOnOff and runs the function train() where the
 * main functionality is implemented.
 *
 * This examples can generate the matrix needed for the density estimation.
 * However the the matrix can be build and decomposed with the class "buildMats.cpp"
 * This allows to save the matrix decomposition and use it multiple times.
 * 
 */



int main() {
  /**
   * Get the training, test and validation data
   */
  std::string filename = "../../datasets/mnist/mnist_[2, 5, 7]_28x28.train.arff";
  // load training samples
  std::cout << "# loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset trainDataset = sgpp::datadriven::ARFFTools::readARFF(filename);

  filename = "../../datasets/mnist/mnist_[2, 5, 7]_28x28.t10k.arff";
  // load test samples
  std::cout << "# loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset testDataset = sgpp::datadriven::ARFFTools::readARFF(filename);

  // if fixed validation data should be used (required for convergence
  // monitor):
  /*filename = "";  // specify file containing validation data here
  // load validation samples
  std::cout << "# loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset validationDataset =
      sgpp::datadriven::ARFFTools::readARFF(filename); */

  /**
   * Specify the number of classes and the corresponding class labels.
   */
  size_t classNum = 3;
  sgpp::base::DataVector classLabels(classNum);
  classLabels[0] = 2;
  classLabels[1] = 5;
  classLabels[2] = 7;
  /*
  classLabels[2] = 2;
  classLabels[3] = 3;
  classLabels[4] = 4;
  classLabels[5] = 5;
  classLabels[6] = 6;
  classLabels[7] = 7;
  classLabels[8] = 8;
  classLabels[9] = 9;
  */
  /**
   * The grid configuration.
   * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   * These Setting get overwritten if a matrixfile is provided
   * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   */
  std::cout << "# create grid config" << std::endl;
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = trainDataset.getDimension();
  gridConfig.level_ = 3;
  // gridConfig.type_ = sgpp::base::GridType::Linear;
  gridConfig.type_ = sgpp::base::GridType::ModLinear;

  /**
   * Configure regularization.
   */
  std::cout << "# create regularization config" << std::endl;
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
  // initial regularization parameter lambda
  regularizationConfig.lambda_ = 0.01;
  /**
   * Select the desired decomposition type for the offline step.
   * Note: Refinement/Coarsening only possible for Cholesky decomposition.
   */
  sgpp::datadriven::MatrixDecompositionType dt;
  std::string decompType;
  // choose "LU decomposition"
  // dt = MatrixDecompositionType::DBMatDecompLU;
  // decompType = "LU decomposition";
  // choose"Eigen decomposition"
  // dt = MatrixDecompositionType::DBMatDecompEigen;
  // decompType = "Eigen decomposition";
  // choose "Cholesky decomposition"
  dt = sgpp::datadriven::MatrixDecompositionType::Chol;
  decompType = "Cholesky decomposition";
  //      dt = sgpp::datadriven::MatrixDecompositionType::IChol;
  //      decompType = "Incomplete Cholesky decomposition";
  //    dt = sgpp::datadriven::MatrixDecompositionType::DenseIchol;
  //    decompType = "Incomplete Cholesky decomposition on Dense Matrix";
  std::cout << "Decomposition type: " << decompType << std::endl;
  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  densityEstimationConfig.decomposition_ = dt;
  densityEstimationConfig.iCholSweepsDecompose_ = 2;
  densityEstimationConfig.iCholSweepsRefine_ = 2;

  /**
   * Configure adaptive refinement (if Cholesky is chosen). As refinement
   * monitor the periodic monitor or the convergence monitor
   * can be chosen. Possible refinement indicators are
   * surplus refinement, data-based refinement, zero-crossings-based
   * refinement.
   */
  std::cout << "# create adaptive refinement configuration" << std::endl;
  std::string refMonitor;
  // select periodic monitor - perform refinements in fixed intervals
  refMonitor = "periodic";
  size_t refPeriod = 1;  // the refinement interval
  // select convergence monitor - perform refinements if algorithm has
  // converged
  // (convergence measured with respect to changes of the classification
  // accuracy)
  // refMonitor = "convergence";
  // the convergence threshold
  double accDeclineThreshold = 0.001;
  // number of accuracy measurements which
  // are considered for convergence check
  size_t accDeclineBufferSize = 140;
  // minimum number of iterations before next refinement
  // is allowed to be performed
  size_t minRefInterval = 10;
  std::cout << "Refinement monitor: " << refMonitor << std::endl;
  std::string refType;
  // select surplus refinement
  // refType = "surplus";
  // select data-based refinement
  // refType = "data";
  // select zero-crossings-based refinement
  refType = "zero";
  std::cout << "Refinement type: " << refType << std::endl;
  sgpp::base::AdpativityConfiguration adaptConfig;
  /**
   * Specify number of refinement steps and the max number
   * of grid points to refine each step.
   */
  adaptConfig.numRefinements_ = 10;
  adaptConfig.noPoints_ = 5;
  adaptConfig.threshold_ = 0.0;  // only required for surplus refinement

  // initial weighting factor
  double beta = 0.0;
  // specify if prior should be used to predict class labels
  bool usePrior = false;


  std::string matrixfile = "mats/28x28_ModLin_NN_Inter_lvl3_Chol.out";

  /**
   * Create the learner.
   */
  std::cout << "# create learner" << std::endl;
  // Use precomputed matrix
  sgpp::datadriven::LearnerSGDEOnOff learner(gridConfig, adaptConfig, regularizationConfig,
                                             densityEstimationConfig, trainDataset, testDataset,
                                             nullptr, classLabels, classNum, usePrior, beta,
                                             matrixfile);


  // Build matrix
  // sgpp::datadriven::LearnerSGDEOnOff learner(gridConfig, adaptConfig, regularizationConfig,
  //    densityEstimationConfig, trainDataset, testDataset, nullptr,
  //                                           classLabels, classNum, usePrior, beta);
  std::cout << "learner set up!" << std::endl;

  /**
   * Configure cross-validation.
   * Set enableCv=true to perform cross-validation
   * during the learning process.
   */
  bool enableCv = false;
  // set cv configuration if cv enabled
  size_t nextCvStep = 50;
  double cvLambdaStart = 1e-1;
  double cvLambdaEnd = 1e-10;
  int cvLambdaSteps = 10;
  bool cvLogScale = true;
  sgpp::base::DataMatrix* cvTestData = &testDataset.getData();
  sgpp::base::DataMatrix* cvTestDataRes = nullptr;  // needed?
  learner.setCrossValidationParameters(cvLambdaSteps, cvLambdaStart, cvLambdaEnd, cvTestData,
                                       cvTestDataRes, cvLogScale);

  // specify batch size
  // (set to 1 for processing only a single data point each iteration)
  size_t batchSize = 17644;
  // specify max number of passes over traininig data set
  size_t maxDataPasses = 1;

  /**
   * Learn the data.
   */
  std::cout << "# start to train the learner" << std::endl;
  learner.train(batchSize, maxDataPasses, refType, refMonitor, refPeriod, accDeclineThreshold,
                accDeclineBufferSize, minRefInterval, enableCv, nextCvStep);

  // learner.storeResults();

  auto xTest = testDataset.getData();
  auto yTest = testDataset.getTargets();

  std::cout << "Start prediction..." << std::endl;

  sgpp::base::DataVector predictTest(xTest.size());
  learner.predict(xTest, predictTest);

  std::cout << "Build confusion matrix \n real\\pred" << std::endl;

  classNum = 10;
  size_t *confusion = new size_t[classNum*classNum];

  for (size_t real = 0; real < classNum; real++) {
    for (size_t pred = 0; pred < classNum; pred++) {
      confusion[pred + real*classNum] = 0;
    }
  }

  for (size_t i = 0; i < testDataset.getNumberInstances(); i++) {
    confusion[(size_t)predictTest.get(i)+classNum*(size_t)yTest.get(i)]++;
  }

  for (size_t real = 0; real < classNum; real++) {
    for (size_t pred = 0; pred < classNum; pred++) {
      std::cout << confusion[pred+real*classNum] << "\t";
    }
    std::cout << std::endl;
  }

  /**
   * Accuracy on test data.
   */
  double acc = learner.getAccuracy();
  std::cout << "# accuracy (test data): " << acc << std::endl;
/*
  double trainerr = learner.getError(trainDataset);
  std::cout << "Train error: " << trainerr << std::endl;
*/
  // store results (classified data, grids, density functions)
  // learner.storeResults();
}

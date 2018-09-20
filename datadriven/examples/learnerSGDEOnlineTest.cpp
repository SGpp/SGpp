// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/application/LearnerSGDE.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>

#include <string>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

/**
 * \page example_learnerSGDEOnlineTest_cpp Learner SGDE Online
 * This example shows how to perform online-classification using sparse
 * grid density estimation and conjugate gradients method. It creates an
 * instance of LearnerSGDE and runs the function trainOnline() where the
 * main functionality is implemented.
 *
 * Currently, only binary classification with class labels -1 and 1 is possible.
 *
 * The example provides the option to execute several runs over differently
 * ordered data and perform a 5-fold cross-validation within each run.
 * Therefore,
 * already randomly ordered and partitioned data is required.
 * Average results from several runs might be more reliable in an
 * online-learning
 * scenario, because the ordering of the data points seen by the learner
 * can affect the result.
 */

int main() {
  /**
   * Specify the number of runs to perform.
   * If only one specific example should be executed, set
   * totalSets=1.
   */
  size_t totalSets = 1;
  size_t totalFolds = 1;  // set to 5 to perform 5-fold cv
  double avgError = 0.0;
  double avgErrorFolds = 0.0;
  for (size_t numSets = 0; numSets < totalSets; numSets++) {
    /**
     * A vector to compute average classification error throughout
     * the learning process. The length of the vector determines
     * the total number of error observations.
     */
    sgpp::base::DataVector avgErrorsFolds(51, 0.0);

    for (size_t numFolds = 0; numFolds < totalFolds; numFolds++) {
      /**
       * Get the training, test and validation data
       */
      std::string filename = "../../datasets/ripley/ripleyGarcke.train.arff";
      // load training samples
      std::cout << "# loading file: " << filename << std::endl;
      sgpp::datadriven::Dataset trainDataset = sgpp::datadriven::ARFFTools::readARFFFromFile(filename);
      sgpp::base::DataMatrix& trainData = trainDataset.getData();
      // extract training classes
      sgpp::base::DataVector& trainLabels = trainDataset.getTargets();

      filename = "../../datasets/ripley/ripleyGarcke.test.arff";
      // load test samples
      std::cout << "# loading file: " << filename << std::endl;
      sgpp::datadriven::Dataset testDataset = sgpp::datadriven::ARFFTools::readARFFFromFile(filename);
      sgpp::base::DataMatrix& testData = testDataset.getData();
      // extract test classes
      sgpp::base::DataVector& testLabels = testDataset.getTargets();

      sgpp::base::DataMatrix* validData = nullptr;
      sgpp::base::DataVector* validLabels = nullptr;
      // if fixed validation data should be used (required for convergence
      // monitor):
      // filename = "";  // specify file containing validation data here
      // load validation samples
      // std::cout << "# loading file: " << filename << std::endl;
      // sgpp::datadriven::Dataset valDataset =
      //    sgpp::datadriven::ARFFTools::readARFF(filename);
      // validData = &valDataset.getData();
      // extract validation classes
      // validLabels = &valDataset.getTargets();

      /**
       * Specify the ocurring class labels.
       */
      size_t classNum = 2;
      sgpp::base::DataVector classLabels(classNum);
      classLabels[0] = -1;
      classLabels[1] = 1;

      /**
       * The grid configuration.
       */
      std::cout << "# create grid config" << std::endl;
      sgpp::base::RegularGridConfiguration gridConfig;
      gridConfig.dim_ = trainDataset.getDimension();
      gridConfig.level_ = 3;
      gridConfig.type_ = sgpp::base::GridType::Linear;
      // gridConfig.type_ = sgpp::base::GridType::ModLinear;

      /**
       * Configure adaptive refinement. As refinement
       * monitor the periodic monitor or the convergence monitor
       * can be chosen. Possible refinement indicators are
       * surplus refinement, data-based refinement, zero-crossings-based
       * refinement.
       */
      std::cout << "# create adaptive refinement config" << std::endl;
      std::string refMonitor;
      // select periodic monitor - perform refinements in fixed intervals
      refMonitor = "periodic";
      size_t refPeriod = 40;  // the refinement interval
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
      /**
       * Specify number of refinement steps and the max number
       * of grid points to refine each step.
       */
      sgpp::base::AdaptivityConfiguration adaptConfig;
      adaptConfig.numRefinements_ = 2;
      adaptConfig.noPoints_ = 7;
      adaptConfig.threshold_ = 0.0;  // only required for surplus refinement

      /**
       * Configure the CG solver. Note that the max number of
       * iterations should be limited in order to obtain
       * feasible runtimes, especially for large grids.
       */
      std::cout << "# create solver config" << std::endl;
      sgpp::solver::SLESolverConfiguration solverConfig;
      solverConfig.type_ = sgpp::solver::SLESolverType::CG;
      solverConfig.maxIterations_ = 20;
      solverConfig.eps_ = 1e-10;
      solverConfig.threshold_ = 1e-10;

      /**
       * Configure regularization.
       */
      std::cout << "# create regularization config" << std::endl;
      sgpp::datadriven::RegularizationConfiguration regularizationConfig;
      regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
      // regularizationConfig.type_ =
      // sgpp::datadriven::RegularizationType::Laplace;

      /**
       * Configure cross-validation.
       */
      std::cout << "# create cross-validation config" << std::endl;
      sgpp::datadriven::CrossvalidationConfiguration crossvalidationConfig;
      crossvalidationConfig.lambda_ = 0.01;
      crossvalidationConfig.enable_ = false;  // set 'true' to perform cv

      crossvalidationConfig.kfold_ = 5;
      crossvalidationConfig.lambdaStart_ = 1e-1;
      crossvalidationConfig.lambdaEnd_ = 1e-10;
      crossvalidationConfig.lambdaSteps_ = 5;
      crossvalidationConfig.logScale_ = true;
      crossvalidationConfig.shuffle_ = true;
      crossvalidationConfig.seed_ = 1234567;
      crossvalidationConfig.silent_ = true;

      /**
       * Create the learner.
       */
      std::cout << "# creating the learner" << std::endl;
      sgpp::datadriven::LearnerSGDE learner(gridConfig, adaptConfig, solverConfig,
                                            regularizationConfig, crossvalidationConfig);
      learner.initialize(trainData);

      // specify if prior should be used to predict class labels
      bool usePrior = false;
      // specify max number of passes over traininig data set
      size_t maxDataPasses = 2;

      /**
       * Learn the data.
       */
      std::cout << "# start to train the learner" << std::endl;
      learner.trainOnline(trainLabels, testData, testLabels, validData, validLabels, classLabels,
                          maxDataPasses, refType, refMonitor, refPeriod, accDeclineThreshold,
                          accDeclineBufferSize, minRefInterval, usePrior);

      std::cout << "# finished training" << std::endl;

      /**
       * Accuracy on test and current training data.
       */
      double accTrain = learner.getAccuracy(trainData, trainLabels, 0.0);
      std::cout << "Acc (train): " << accTrain << std::endl;
      double accTest = learner.getAccuracy(testData, testLabels, 0.0);
      std::cout << "Acc (test): " << accTest << std::endl;

      // store results (classified data, grids, density functions)
      // learner.storeResults(testData);

      avgErrorFolds += learner.error;
      avgErrorsFolds.add(learner.avgErrors);
    }
    avgErrorFolds = avgErrorFolds / static_cast<double>(totalFolds);
    /**
     * Average accuracy on test data reagarding 5-fold cv.
     */
    if ((totalSets > 1) && (totalFolds > 1)) {
      std::cout << "Average accuracy on test data (set " + std::to_string(numSets + 1) + "): "
                << (1.0 - avgErrorFolds) << "\n";
    }
    ///
    avgError += avgErrorFolds;
    avgErrorFolds = 0.0;

    avgErrorsFolds.mult(1.0 / static_cast<double>(totalFolds));

    // write error evaluation to csv file
    // std::ofstream output;
    // output.open("SGDE_avg_classification_error_"+std::to_string(numSets+1)+".csv");
    // if (output.fail()) {
    //  std::cout << "failed to create csv file!" << std::endl;
    // }
    // else {
    //  for (size_t i = 0; i < avgErrorsFolds.getSize(); i++) {
    //    output << avgErrorsFolds.get(i) << ";" << std::endl;
    // }
    //  output.close();
    // }
  }
  ///
}

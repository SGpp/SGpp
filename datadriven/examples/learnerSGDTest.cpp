// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/application/LearnerSGD.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/globaldef.hpp>
#include <string>

/**
 * \page example_learnerSGDTest_cpp Learner SGD
 * This example shows how to perform online-classification using sparse
 * grids and averaged stochastic gradient descent method. It creates an
 * instance of LearnerSGD and runs the function train() where the
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
      std::string filename = "../datasets/ripley/ripleyGarcke.train.arff";
      // load training samples
      std::cout << "# loading file: " << filename << std::endl;
      sgpp::datadriven::Dataset trainDataset =
          sgpp::datadriven::ARFFTools::readARFFFromFile(filename);
      sgpp::base::DataMatrix& trainData = trainDataset.getData();
      // extract train classes
      sgpp::base::DataVector& trainLabels = trainDataset.getTargets();

      filename = "../datasets/ripley/ripleyGarcke.test.arff";
      // load test samples
      std::cout << "# loading file: " << filename << std::endl;
      sgpp::datadriven::Dataset testDataset =
          sgpp::datadriven::ARFFTools::readARFFFromFile(filename);
      sgpp::base::DataMatrix& testData = testDataset.getData();
      // extract test classes
      sgpp::base::DataVector& testLabels = testDataset.getTargets();

      sgpp::base::DataMatrix* validData = nullptr;
      sgpp::base::DataVector* validLabels = nullptr;
      bool useValidData = false;
      // if fixed validation data should be used (required for convergence
      // monitor):
      // filename = "";  // specify file containing validation data here
      // load validation samples
      // std::cout << "# loading file: " << filename << std::endl;
      // sgpp::datadriven::Dataset valDataset =
      //    sgpp::datadriven::ARFFTools::readARFF(filename);
      // validData = &(valDataset.getData());
      // extract validation classes
      // validLabels = &(valDataset.getTargets());
      // useValidData = true;

      /**
       * The grid configuration.
       */
      std::cout << "# creating grid config" << std::endl;
      sgpp::base::RegularGridConfiguration gridConfig;
      gridConfig.dim_ = trainDataset.getDimension();
      gridConfig.level_ = 3;
      // gridConfig.type_ = sgpp::base::GridType::Linear;
      gridConfig.type_ = sgpp::base::GridType::ModLinear;

      /**
       * Configure adaptive refinement. As refinement
       * monitor the periodic monitor or the convergence monitor
       * can be chosen. Possible refinement indicators are
       * predictive refinement and impurity-based refinement.
       */
      std::cout << "# create adaptive refinement config" << std::endl;
      std::string refMonitor;
      // select periodic monitor - perform refinements in fixed intervals
      refMonitor = "periodic";
      size_t refPeriod = 40;  // the refinement interval
      // select convergence monitor - perform refinements if algorithm has
      // converged
      // (convergence measured with respect to MSE or Hinge loss observations)
      // refMonitor = "convergence";
      // the convergence threshold
      double errorDeclineThreshold = 0.0005;
      // number of error measurements which
      // are considered for convergence check
      size_t errorDeclineBufferSize = 100;
      // minimum number of iterations before next refinement
      // is allowed to be performed
      size_t minRefInterval = 0;
      std::cout << "Refinement monitor: " << refMonitor << std::endl;
      std::string refType;
      // select predictive refinement
      // refType = "predictive";
      // select impurity-based refinement
      refType = "impurity";
      std::cout << "Refinement type: " << refType << std::endl;
      sgpp::base::AdaptivityConfiguration adaptivityConfig;
      /**
       * Specify number of refinement steps and the max number
       * of grid points to refine each step.
       */
      adaptivityConfig.numRefinements_ = 3;
      adaptivityConfig.numRefinementPoints_ = 5;
      adaptivityConfig.refinementThreshold_ = 0.0;

      // additional parameters:

      // specify max number of passes over traininig data set
      size_t maxDataPasses = 2;
      // regularization parameter
      double lambda = 0.01;
      // initial learning rate
      double gamma = 0.25;

      // specify number of data points to compute error contributions
      // for predictive refinement (if no validation data provided)
      size_t batchSize = 50;
      if (validData != nullptr) {
        batchSize = validData->getNrows();
      }

      /**
       * Create the learner.
       */
      std::cout << "# creating the learner" << std::endl;
      sgpp::datadriven::LearnerSGD learner(gridConfig, adaptivityConfig, trainData, trainLabels,
                                           testData, testLabels, validData, validLabels, lambda,
                                           gamma, batchSize, useValidData);

      // initialize learner (create grid etc.)
      learner.initialize();

      /**
       * Learn the data.
       */
      std::cout << "# start to train the learner" << std::endl;
      learner.train(maxDataPasses, refType, refMonitor, refPeriod, errorDeclineThreshold,
                    errorDeclineBufferSize, minRefInterval);

      // store results (classified data, grid, function evaluations)
      // learner.storeResults(testData);

      /**
       * Accuracy on test and current training data.
       */
      double accTrain = learner.getAccuracy(trainData, trainLabels, 0.0);
      std::cout << "Acc (train): " << accTrain << std::endl;
      double accTest = learner.getAccuracy(testData, testLabels, 0.0);
      std::cout << "Acc (test): " << accTest << std::endl;

      avgErrorFolds += learner.error;
      avgErrorsFolds.add(learner.avgErrors);
    }
    avgErrorFolds = avgErrorFolds / static_cast<double>(totalFolds);
    if ((totalSets > 1) && (totalFolds > 1)) {
      /**
       * Average accuracy on test data reagarding 5-fold cv.
       */
      std::cout << "Average accuracy on test data (set " + std::to_string(numSets + 1) + "): "
                << (1.0 - avgErrorFolds) << std::endl;
    }
    avgError += avgErrorFolds;
    avgErrorFolds = 0.0;

    avgErrorsFolds.mult(1.0 / static_cast<double>(totalFolds));

    // write error evaluation to csv file
    // std::ofstream output;
    // output.open("ASGD_avg_classification_error_"+std::to_string(numSets+1)+".csv");
    // if (output.fail()) {
    //  std::cout << "failed to create .csv file!" << std::endl;
    //}
    // else {
    //  for (size_t i = 0; i < avgErrorsFolds.getSize(); i++) {
    //    output << avgErrorsFolds.get(i) << ";" << std::endl;
    //  }
    //  output.close();
    //}
  }
  // avgError = avgError / static_cast<double>(totalSets);
  // std::cout << "Average accuracy on test data: " << (1.0 - avgError) <<
  // std::endl;
  ///
}

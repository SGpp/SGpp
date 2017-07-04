// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp>
#include <sgpp/datadriven/application/LearnerSGDEOnOff.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>

#include <string>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

/**
 * This example shows how to perform offline/online-classification using sparse
 * grid density estimation and matrix decomposition methods. It creates an
 * instance of LearnerSGDEOnOff and runs the function train() where the
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
      sgpp::datadriven::Dataset trainDataset = sgpp::datadriven::ARFFTools::readARFF(filename);

      filename = "../../datasets/ripley/ripleyGarcke.test.arff";
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
       * Configure regularization.
       */
      std::cout << "# create regularization config" << std::endl;
      sgpp::datadriven::RegularizationConfiguration regularizationConfig;
      regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Identity;

      /**
       * Select the desired decomposition type for the offline step.
       * Note: Refinement/Coarsening only possible for Cholesky decomposition.
       */
      sgpp::datadriven::DBMatDecompostionType dt;
      std::string decompType;
      // choose "LU decomposition"
      // dt = DBMatDecompostionType::DBMatDecompLU;
      // decompType = "LU decomposition";
      // choose"Eigen decomposition"
      // dt = DBMatDecompostionType::DBMatDecompEigen;
      // decompType = "Eigen decomposition";
      // choose "Cholesky decomposition"
      //      dt = sgpp::datadriven::DBMatDecompostionType::Chol;
      //      decompType = "Cholesky decomposition";
      //      dt = sgpp::datadriven::DBMatDecompostionType::IChol;
      //      decompType = "Incomplete Cholesky decomposition";
      dt = sgpp::datadriven::DBMatDecompostionType::DenseIchol;
      decompType = "Incomplete Cholesky decomposition on Dense Matrix";
      std::cout << "Decomposition type: " << decompType << std::endl;

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
      sgpp::base::AdpativityConfiguration adaptConfig;
      /**
       * Specify number of refinement steps and the max number
       * of grid points to refine each step.
       */
      adaptConfig.numRefinements_ = 2;
      adaptConfig.noPoints_ = 7;
      adaptConfig.threshold_ = 0.0;  // only required for surplus refinement

      // initial regularization parameter lambda
      double lambda = 0.01;
      // initial weighting factor
      double beta = 0.0;
      // configuration
      sgpp::datadriven::DBMatDensityConfiguration dconf(gridConfig, adaptConfig,
                                                        regularizationConfig.regType_, lambda, dt);
      // specify if prior should be used to predict class labels
      bool usePrior = false;

      dconf.icholParameters.sweepsDecompose = 2;
      dconf.icholParameters.sweepsRefine = 2;

      /**
       * Create the learner.
       */
      std::cout << "# create learner" << std::endl;
      sgpp::datadriven::LearnerSGDEOnOff learner(dconf, trainDataset, testDataset, nullptr,
                                                 classLabels, classNum, usePrior, beta, lambda);

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
      size_t batchSize = 1;
      // specify max number of passes over traininig data set
      size_t maxDataPasses = 2;

      /**
       * Learn the data.
       */
      std::cout << "# start to train the learner" << std::endl;
      learner.train(batchSize, maxDataPasses, refType, refMonitor, refPeriod, accDeclineThreshold,
                    accDeclineBufferSize, minRefInterval, enableCv, nextCvStep);

      /**
       * Accuracy on test data.
       */
      double acc = learner.getAccuracy();
      std::cout << "# accuracy (test data): " << acc << std::endl;

      // store results (classified data, grids, density functions)
      // learner.storeResults();

      sgpp::base::DataVector tmp{};
      avgErrorFolds += 1.0 - learner.getAccuracy();
      learner.getAvgErrors(tmp);
      avgErrorsFolds.add(tmp);
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

    // write error evaluation to csv-file
    /*std::ofstream output;
    output.open("SGDEOnOff_avg_classification_error_"+std::to_string(numSets+1)+".csv");
    if (output.fail()) {
      std::cout << "failed to create csv file!" << std::endl;
    }
    else {
      for (size_t i = 0; i < avgErrorsFolds.getSize(); i++) {
        output << avgErrorsFolds.get(i) << ";" << std::endl;
      }
      output.close();
    }*/
  }
}

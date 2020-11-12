// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#include <sgpp/datadriven/configuration/DensityEstimationConfiguration.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>

#endif /* USE_GSL */

#ifdef USE_MPI
#include <omp.h>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/LearnerSGDEOnOffParallel.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/MPIMethods.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/RoundRobinScheduler.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <string>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

sgpp::datadriven::Dataset loadDataset(const std::string &filename);
void parseInputValue(char *inputString, size_t &outputValue);
#endif /* USE_MPI */

/**
 * \page example_learnerSGDEOnOffTest_cpp Learner SGDE OnOff
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

int main(int argc, char *argv[]) {
#ifdef USE_MPI
#ifdef USE_GSL

  omp_set_num_threads(1);

  std::cout << "LearnerSGDEOnOffParallelTest" << std::endl;

  if (argc != 5) {
    std::cout << "Usage:" << std::endl
              << "learnerSGDEOnOffParallelTest <trainDataFile> "
              << "<testDataFile> <batchSize> <refPeriod>" << std::endl;
    return -1;
  }
  /**
   * Specify the number of runs to perform.
   * If only one specific example should be executed, set
   * totalSets=1.
   */
  //  size_t totalSets = 1;
  //  size_t totalFolds = 1;  // set to 5 to perform 5-fold cv
  //  double avgError = 0.0;
  //  double avgErrorFolds = 0.0;
  //  for (size_t numSets = 0; numSets < totalSets; numSets++) {
  //    /**
  //     * A vector to compute average classification error throughout
  //     * the learning process. The length of the vector determines
  //     * the total number of error observations.
  //     */
  //    sgpp::base::DataVector avgErrorsFolds(51, 0.0);
  //
  //    for (size_t numFolds = 0; numFolds < totalFolds; numFolds++) {
  /**
   * Get the training, test and validation data
   */
  sgpp::datadriven::Dataset trainDataset = loadDataset(argv[1]);
  sgpp::datadriven::Dataset testDataset = loadDataset(argv[2]);

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
  sgpp::datadriven::RegularizationConfiguration regularizationConfig{};
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
  //      dt = sgpp::datadriven::MatrixDecompositionType::Chol;
  //      decompType = "Cholesky decomposition";
  //      dt = sgpp::datadriven::MatrixDecompositionType::IChol;
  //      decompType = "Incomplete Cholesky decomposition";
  dt = sgpp::datadriven::MatrixDecompositionType::DenseIchol;
  decompType = "Incomplete Cholesky decomposition on Dense Matrix";
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

  size_t refPeriod = 0;  // the refinement interval
  parseInputValue(argv[4], refPeriod);
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
  sgpp::base::AdaptivityConfiguration adaptivityConfig;
  /**
   * Specify number of refinement steps and the max number
   * of grid points to refine each step.
   */
  adaptivityConfig.numRefinements_ = 2;
  adaptivityConfig.numRefinementPoints_ = 7;
  adaptivityConfig.refinementThreshold_ = 0.0;  // only required for surplus refinement

  // initial weighting factor
  double beta = 0.0;
  // specify if prior should be used to predict class labels
  bool usePrior = false;

  // specify batch size
  // (set to 1 for processing only a single data point each iteration)
  size_t batchSize = 0;
  parseInputValue(argv[3], batchSize);

  // Create the MPI Task Scheduling using the round robin algorithm
  sgpp::datadriven::RoundRobinScheduler scheduler(batchSize);

  /**
   * Create the learner.
   */
  std::cout << "# create learner" << std::endl;
  sgpp::datadriven::LearnerSGDEOnOffParallel learner(
      gridConfig, adaptivityConfig, regularizationConfig, densityEstimationConfig, trainDataset,
      testDataset, nullptr, classLabels, classNum, usePrior, beta, scheduler);

  // specify max number of passes over traininig data set
  size_t maxDataPasses = 1;

  /**
   * Learn the data.
   */

  MPI_Barrier(MPI_COMM_WORLD);

  std::cout << "# start to train the learner" << std::endl;
  sgpp::base::SGppStopwatch stopwatch;

  //    CALLGRIND_START_INSTRUMENTATION;

  stopwatch.start();
  learner.trainParallel(batchSize, maxDataPasses, refType, refMonitor, refPeriod,
                        accDeclineThreshold, accDeclineBufferSize, minRefInterval);
  double deltaTime = stopwatch.stop();

  //    CALLGRIND_STOP_INSTRUMENTATION;

  MPI_Barrier(MPI_COMM_WORLD);

  /**
   * Accuracy on test data.
   */
  double acc = learner.getAccuracy();
  if (sgpp::datadriven::MPIMethods::isMaster()) {
    std::cout << "# accuracy (test data): " << acc << std::endl;
    std::cout << "# delta time training: " << deltaTime << std::endl;

  } else {
    std::cout << "# accuracy (client, test data): " << acc << std::endl;
    std::cout << "# delta time training (client): " << deltaTime << std::endl;
  }
// store results (classified data, grids, density functions)
// learner.storeResults();

//        sgpp::base::DataVector tmp;
//        avgErrorFolds += 1.0 - learner.getAccuracy();
//        learner.getAvgErrors(tmp);
//        avgErrorsFolds.add(tmp);
//    }
//    avgErrorFolds = avgErrorFolds / static_cast<double>(totalFolds);
//    if ((totalSets > 1) && (totalFolds > 1)) {
//      /**
//       * Average accuracy on test data reagarding 5-fold cv.
//       */
//      std::cout << "Average accuracy on test data (set " + std::to_string(numSets + 1) + "): "
//                << (1.0 - avgErrorFolds) << std::endl;
//    }
//    avgError += avgErrorFolds;
//    avgErrorFolds = 0.0;
//    avgErrorsFolds.mult(1.0 / static_cast<double>(totalFolds));

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
//  }
#else
  std::cout << "GSL not enabled at compile time" << std::endl;
#endif  // USE_GSL
#endif  /* USE_MPI */
  ///
}

#ifdef USE_MPI
sgpp::datadriven::Dataset loadDataset(const std::string &filename) {
  // load test samples
  std::cout << "# loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset dataset = sgpp::datadriven::ARFFTools::readARFFFromFile(filename);

  if (dataset.getDimension() <= 0) {
    std::cout << "# Failed to read dataset! " << filename << std::endl;
    exit(-1);
  } else {
    std::cout << "# dataset dimensionality: " << dataset.getDimension() << std::endl;
  }
  return dataset;
}

void parseInputValue(char *inputString, size_t &outputValue) {
  std::basic_stringstream<char> argumentParser = std::stringstream(inputString);
  argumentParser >> outputValue;
  if (argumentParser.fail()) {
    throw sgpp::base::application_exception("Failed to parse parameter");
  }
}
#endif /* USE_MPI */

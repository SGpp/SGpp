/*
 * learnerOnOff.cpp
 *
 *  Created on: Apr 27, 2017
 *      Author: milettri
 */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#ifdef USE_GSL
#include <sgpp/datadriven/application/LearnerSGDEOnOff.hpp>
#include <sgpp/datadriven/configuration/DensityEstimationConfiguration.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#endif /* USE_GSL */
#include <sgpp/datadriven/tools/ARFFTools.hpp>

#include <string>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

int main() {
#ifdef USE_GSL

  /**
   * Get the training, test and validation data
   */
  auto filename = "../../datasets/ripley/ripleyGarcke.train.arff";
  // auto filename = "dr10_train.arff";
  // load training samples
  std::cout << "# loading file: " << filename << std::endl;
  auto trainDataset = sgpp::datadriven::ARFFTools::readARFF(filename);
  filename = "../../datasets/ripley/ripleyGarcke.test.arff";
  // load test samples
  // filename = "dr10_test.arff";
  std::cout << "# loading file: " << filename << std::endl;
  auto testDataset = sgpp::datadriven::ARFFTools::readARFF(filename);

  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.type_ = sgpp::base::GridType::Linear;
  gridConfig.dim_ = trainDataset.getDimension();
  gridConfig.level_ = 2;

  sgpp::base::AdpativityConfiguration adaptConfig;
  adaptConfig.numRefinements_ = 10;
  adaptConfig.noPoints_ = 1;
  adaptConfig.threshold_ = 0.0;

  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;
  regularizationConfig.lambda_ = 5 * 10e-5;

  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  densityEstimationConfig.iCholSweepsDecompose_ = 2;
  densityEstimationConfig.iCholSweepsRefine_ = 2;
  densityEstimationConfig.iCholSweepsSolver_ = 2;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::DenseIchol;
  auto decompType = "Incomplete Cholesky decomposition on Dense Matrix";
  //  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::Chol;
  //  auto decompType = "Cholesky decomposition";
  std::cout << "Decomposition type: " << decompType << std::endl;

  densityEstimationConfig.iCholSweepsDecompose_ = 8;
  densityEstimationConfig.iCholSweepsRefine_ = 8;

  /**
   * Specify the number of classes and the corresponding class labels.
   */
  size_t classNum = 2;
  sgpp::base::DataVector classLabels(classNum);
  classLabels[0] = -1;
  classLabels[1] = 1;

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
  size_t minRefInterval = 1;
  std::cout << "Refinement monitor: " << refMonitor << std::endl;
  std::string refType;
  // select surplus refinement
  refType = "surplus";
  // select data-based refinement
  // refType = "data";
  // select zero-crossings-based refinement
  //  refType = "zero";
  std::cout << "Refinement type: " << refType << std::endl;

  // initial weighting factor
  double beta = 0.0;
  // configuration

  // specify if prior should be used to predict class labels
  bool usePrior = false;

  /**
   * Create the learner.
   */
  std::cout << "# create learner" << std::endl;
  sgpp::datadriven::LearnerSGDEOnOff learner(gridConfig, adaptConfig, regularizationConfig,
                                             densityEstimationConfig, trainDataset, testDataset,
                                             nullptr, classLabels, classNum, usePrior, beta);

  /**
   * Learn the data.
   */
  std::cout << "# start to train the learner" << std::endl;
  learner.train(trainDataset.getNumberInstances(), 10, refType, refMonitor, refPeriod,
                accDeclineThreshold, accDeclineBufferSize, minRefInterval, false, 0);

  /**
   * Accuracy on test data.
   */
  double acc = learner.getAccuracy();
  std::cout << "# accuracy (test data): " << acc << std::endl;
  std::cout << "#####################################################################\n";

// store results (classified data, grids, density functions)
// learner.storeResults();

#endif /* USE_GSL */
}

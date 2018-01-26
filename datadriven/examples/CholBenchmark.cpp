/*
 * CholBenchmark.cpp
 *
 *  Created on: Apr 25, 2017
 *      Author: milettri
 */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#ifdef USE_GSL
#include <sgpp/datadriven/algorithm/DBMatOfflineDenseIChol.hpp>
#include <sgpp/datadriven/application/LearnerSGDEOnOff.hpp>
#include <sgpp/datadriven/configuration/DecompositionConfiguration.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>

#include <chrono>
#include <string>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

#endif /* USE_GSL */

int main() {
#ifdef USE_GSL

  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = 4;
  gridConfig.level_ = 5;

  sgpp::base::AdpativityConfiguration adaptConfig;

  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.lambda_ = 0;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;

  sgpp::datadriven::DecompositionConfiguration decompositionConfig;
  decompositionConfig.iCholSweepsDecompose_ = 2;
  decompositionConfig.type_ = sgpp::datadriven::DBMatDecompostionType::Chol;

  auto decompType = "Incomplete Cholesky decomposition on Dense Matrix";
  std::cout << "Decomposition type: " << decompType << std::endl;

  sgpp::datadriven::DBMatOfflineDenseIChol offline(gridConfig, adaptConfig,
                                                   regularizationConfig, decompositionConfig);
  // sgpp::datadriven::DBMatOfflineChol offline(gridConfig, adaptConfig,
  //           regularizationConfig, decompositionConfig);

  offline.buildMatrix();

  auto begin = std::chrono::high_resolution_clock::now();
  offline.decomposeMatrix();
  auto end = std::chrono::high_resolution_clock::now();
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms"
            << std::endl;
#endif /* USE_GSL */
}

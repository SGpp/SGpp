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
#include <sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineDenseIChol.hpp>
#include <sgpp/datadriven/application/LearnerSGDEOnOff.hpp>

#include <chrono>
#include <string>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

#endif /* USE_GSL */

int main() {
#ifdef USE_GSL

  sgpp::datadriven::DBMatDensityConfiguration config;
  config.grid_dim_ = 4;
  config.grid_level_ = 5;
  config.lambda_ = 0;
  config.regularization_ = sgpp::datadriven::RegularizationType::Identity;
  config.icholParameters.sweepsDecompose = 2;

  config.decomp_type_ = sgpp::datadriven::DBMatDecompostionType::Chol;
  auto decompType = "Incomplete Cholesky decomposition on Dense Matrix";
  std::cout << "Decomposition type: " << decompType << std::endl;

  sgpp::datadriven::DBMatOfflineDenseIChol offline(config);
  // sgpp::datadriven::DBMatOfflineChol offline(config);

  offline.buildMatrix();

  auto begin = std::chrono::high_resolution_clock::now();
  offline.decomposeMatrix();
  auto end = std::chrono::high_resolution_clock::now();
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms"
            << std::endl;
#endif /* USE_GSL */
}

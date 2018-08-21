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
#include <sgpp/datadriven/configuration/DensityEstimationConfiguration.hpp>
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

  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  densityEstimationConfig.iCholSweepsDecompose_ = 2;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::Chol;

  auto decompType = "Incomplete Cholesky decomposition on Dense Matrix";
  std::cout << "Decomposition type: " << decompType << std::endl;

  std::unique_ptr<sgpp::base::Grid> grid;
  if (gridConfig.type_ == sgpp::base::GridType::ModLinear) {
    grid =
        std::unique_ptr<sgpp::base::Grid>{sgpp::base::Grid::createModLinearGrid(gridConfig.dim_)};
  } else if (gridConfig.type_ == sgpp::base::GridType::Linear) {
    grid = std::unique_ptr<sgpp::base::Grid>{sgpp::base::Grid::createLinearGrid(gridConfig.dim_)};
  } else {
    return 1;
  }
  sgpp::datadriven::DBMatOfflineDenseIChol offline;
  // sgpp::datadriven::DBMatOfflineChol offline(gridConfig, adaptConfig,
  //           regularizationConfig, densityEstimationConfig);

  offline.buildMatrix(grid.get(), regularizationConfig);
  auto begin = std::chrono::high_resolution_clock::now();
  offline.decomposeMatrix(regularizationConfig, densityEstimationConfig);
  auto end = std::chrono::high_resolution_clock::now();
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms"
            << std::endl;
#endif /* USE_GSL */
}

/*
 * iCholConvergece.cpp
 *
 *  Created on: Apr 25, 2017
 *      Author: milettri
 */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifdef USE_GSL

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/algorithm/DBMatOfflineDenseIChol.hpp>
#include <sgpp/datadriven/application/LearnerSGDEOnOff.hpp>
#include <sgpp/datadriven/configuration/DensityEstimationConfiguration.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/algorithm/GridFactory.hpp>

#include <omp.h>

#include <chrono>
#include <string>
#include <vector>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
#endif /* USE_GSL */

int main() {
#ifdef USE_GSL
  // ###################################################################################
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = 4;
  gridConfig.level_ = 5;

  sgpp::base::AdpativityConfiguration adaptConfig;

  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.lambda_ = 10e-4;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;

  sgpp::datadriven::DensityEstimationConfiguration fullDensityEstimationConfig;
  fullDensityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::Chol;

  sgpp::datadriven::GridFactory gridFactory;
  std::unique_ptr<sgpp::base::Grid> grid = std::unique_ptr<sgpp::base::Grid>{
    gridFactory.createGrid(gridConfig, std::vector<std::vector <size_t>>())
  };

  sgpp::datadriven::DBMatOfflineChol fullOffline;
  fullOffline.buildMatrix(grid.get(), regularizationConfig);
  fullOffline.decomposeMatrix(regularizationConfig, fullDensityEstimationConfig);
  auto& fullMat = fullOffline.getDecomposedMatrix();

  for (auto i = 0u; i < fullMat.getNrows(); i++) {
    for (auto j = i + 1; j < fullMat.getNcols(); j++) {
      fullMat.set(i, j, 0.0);
    }
  }
  // ###################################################################################

  sgpp::datadriven::DensityEstimationConfiguration exactIConfig = fullDensityEstimationConfig;
  fullDensityEstimationConfig.decomposition_ =
  sgpp::datadriven::MatrixDecompositionType::DenseIchol;
  exactIConfig.iCholSweepsDecompose_ = 1;

  sgpp::datadriven::DBMatOfflineDenseIChol exactIOffline;
  exactIOffline.buildMatrix(grid.get(), regularizationConfig);

  auto numThreads = 0;

#pragma omp parallel
  {
#pragma omp single
    { numThreads = omp_get_num_threads(); }
  }
  omp_set_num_threads(1);

  exactIOffline.decomposeMatrix(regularizationConfig, exactIConfig);

  omp_set_num_threads(numThreads);

  auto& exactIMat = exactIOffline.getDecomposedMatrix();

  for (auto i = 0u; i < exactIMat.getNrows(); i++) {
    for (auto j = i + 1; j < exactIMat.getNcols(); j++) {
      exactIMat.set(i, j, 0.0);
    }
  }

  // ###################################################################################
  for (auto i = 1u; i < 10; i++) {
    sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig = exactIConfig;
    densityEstimationConfig.iCholSweepsDecompose_ = i;

    sgpp::datadriven::DBMatOfflineDenseIChol offline;
    offline.buildMatrix(grid.get(), regularizationConfig);
    offline.decomposeMatrix(regularizationConfig, densityEstimationConfig);

    auto& iMat = offline.getDecomposedMatrix();

    for (auto i = 0u; i < iMat.getNrows(); i++) {
      for (auto j = i + 1; j < iMat.getNcols(); j++) {
        iMat.set(i, j, 0.0);
      }
    }

    auto tmpIMat(iMat);
    tmpIMat.sub(exactIMat);
    tmpIMat.abs();
    tmpIMat.sqr();
    std::cout << " ||iichol - ichol|| with " << densityEstimationConfig.iCholSweepsDecompose_
              << " sweeps is: " << std::scientific << std::setprecision(10) << sqrt(tmpIMat.sum())
              << "\n";

    auto tmpIMat2(iMat);
    tmpIMat2.sub(fullMat);
    tmpIMat2.abs();
    tmpIMat2.sqr();
    std::cout << " ||iichol - chol|| with " << densityEstimationConfig.iCholSweepsDecompose_
              << " sweeps is: " << std::scientific << std::setprecision(10) << sqrt(tmpIMat2.sum())
              << "\n";
  }
#endif /* USE_GSL */
}

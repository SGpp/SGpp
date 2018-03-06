// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DBMatDMSOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEOrthoAdapt.hpp>
#include <sgpp/datadriven/configuration/DensityEstimationConfiguration.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>

#include <chrono>
#include <iostream>
#include <vector>

int main() {
  std::cout << "ortho_adapt algorithm benchmarks: \n";

  // config
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = 4;
  gridConfig.level_ = 6;

  sgpp::base::AdpativityConfiguration adaptConfig;

  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.lambda_ = 0.0001;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Identity;

  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::OrthoAdapt;

  size_t number_points_to_refine = 1;
  size_t number_points_to_coarsen = 1;

  std::cout << "dim = " << gridConfig.dim_ << "\n";
  std::cout << "lvl = " << gridConfig.level_ << "\n";
  std::cout << "lambda = " << regularizationConfig.lambda_ << "\n\n";

  // offline phase
  sgpp::datadriven::DBMatOfflineOrthoAdapt offline(gridConfig, adaptConfig,
                                                   regularizationConfig, densityEstimationConfig);

  offline.buildMatrix();
  std::cout << "initial matrix size = " << offline.getDimA();

  std::cout << "\n\ndecomposition took ";
  auto begin = std::chrono::high_resolution_clock::now();
  offline.decomposeMatrix();
  auto end = std::chrono::high_resolution_clock::now();
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms"
            << std::endl;

  // online phase
  sgpp::datadriven::DBMatOnlineDEOrthoAdapt online(offline);

  size_t refine_size = offline.getDimA() + number_points_to_refine;
  size_t coarsen_size = refine_size - number_points_to_coarsen;

  // create random points to refine
  for (size_t i = 0; i < number_points_to_refine; i++) {
    sgpp::base::DataVector vec(refine_size);
    for (size_t j = 0; j < refine_size; j++) {
      double value = (static_cast<double>(rand()) / (RAND_MAX));  // values in [0, 1]
      vec.set(j, value);
    }
    online.add_new_refine_point(vec);
  }

  std::cout << "\nrefining " << number_points_to_refine << " points took ";
  begin = std::chrono::high_resolution_clock::now();
  online.sherman_morrison_adapt(number_points_to_refine, true);
  end = std::chrono::high_resolution_clock::now();
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms"
            << std::endl;

  // create random indices for coarsening
  std::vector<size_t> coarsen_indices = {};
  for (size_t i = 0; i < number_points_to_coarsen; i++) {
    size_t index = ((size_t)rand() % (refine_size - coarsen_size)) + offline.getDimA();
    coarsen_indices.push_back(index);
  }
  // DEBUG:
  // std::cout << "indices to coarsen are: \n";
  // for (auto& i : coarsen_indices) {
  //   std::cout << i << " ";
  // }

  std::cout << "\ncoarsening " << number_points_to_coarsen << " points took ";
  online.sherman_morrison_adapt(0, false, coarsen_indices);
  end = std::chrono::high_resolution_clock::now();
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms"
            << std::endl;

  return 0;
}

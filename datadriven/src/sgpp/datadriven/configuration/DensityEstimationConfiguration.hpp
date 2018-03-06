/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DensityEstimationConfiguration.hpp
 *
 *  Created on: Jan 02, 2018
 *      Author: Kilian Röhner
 */

#pragma once

#include <sgpp/globaldef.hpp>

/**
 * Struct that stores all the configuration information
 * for an offline object for matrix based density estimation
 */

namespace sgpp {
namespace datadriven {

enum class DensityEstimationType { CG, Decomposition };

enum class MatrixDecompositionType { LU, Eigen, Chol, DenseIchol, OrthoAdapt };

struct DensityEstimationConfiguration {
  DensityEstimationType type_;  // Type of density estimation
  MatrixDecompositionType decomposition_;  // Type of matrix decomposition

  // Incomplete Cholesky Decomposition parameters
  size_t iCholSweepsDecompose_ = 4;
  size_t iCholSweepsRefine_ = 4;
  size_t iCholSweepsUpdateLambda_ = 2;
  size_t iCholSweepsSolver_ = 2;
};

}  // namespace datadriven
}  // namespace sgpp

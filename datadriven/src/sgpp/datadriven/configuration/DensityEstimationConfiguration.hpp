// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

/**
 * Struct that stores all the configuration information for an offline object for matrix based
 * density estimation
 */

namespace sgpp {
namespace datadriven {

enum class DensityEstimationType { CG, Decomposition };

enum class MatrixDecompositionType { LU, Eigen, Chol, DenseIchol, OrthoAdapt, SMW_ortho, SMW_chol };

struct DensityEstimationConfiguration {
  // Type of density estimation
  DensityEstimationType type_ = DensityEstimationType::Decomposition;
  // Type of matrix decomposition
  MatrixDecompositionType decomposition_ = MatrixDecompositionType::OrthoAdapt;
  /**
   * Defines whether offline permutation should be used if decomposition allows it.
   */
  bool useOfflinePermutation_ = true;

  // flag for normalization in DBMatOnlineDE
  bool normalize_ = false;

  // Incomplete Cholesky Decomposition parameters
  size_t iCholSweepsDecompose_ = 4;
  size_t iCholSweepsRefine_ = 4;
  size_t iCholSweepsUpdateLambda_ = 2;
  size_t iCholSweepsSolver_ = 2;

  // Partial derivative direction for derivative-based density estimation
  size_t derivDim_ = 0;
  // Relative factor for relative density ratio estimation; value between 0 and 1:
  //    omega = 0 means actual density ratio
  //    omega = 1 means fully damped density (uniform distribution)
  double omega_ = 0;
};

}  // namespace datadriven
}  // namespace sgpp

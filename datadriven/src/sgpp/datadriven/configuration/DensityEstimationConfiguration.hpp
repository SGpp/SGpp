// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

/**
 * Struct that stores all the configuration information
 * for an offline object for matrix based density estimation
 */

namespace sgpp {
namespace datadriven {

enum class DensityEstimationType { CG, Decomposition };

enum class MatrixDecompositionType {
  LU,
  Eigen,
  Chol,
  DenseIchol,
  OrthoAdapt,
  SMW_ortho,
  SMW_chol
};

struct DensityEstimationConfiguration {
  // Type of density estimation
  DensityEstimationType type_ = DensityEstimationType::Decomposition;
  // Type of matrix decomposition
  MatrixDecompositionType decomposition_ = MatrixDecompositionType::OrthoAdapt;
  /**
   * Defines whether offline permutation should be used if decomposition allows
   * it.
   */
  bool useOfflinePermutation = true;

  // flag for normalization in DBMatOnlineDE
  bool normalize_ = false;

  // Incomplete Cholesky Decomposition parameters
  size_t iCholSweepsDecompose_ = 4;
  size_t iCholSweepsRefine_ = 4;
  size_t iCholSweepsUpdateLambda_ = 2;
  size_t iCholSweepsSolver_ = 2;

  /*
  // Debug method to neatly print internal data
  void dumpToStream(std::ostream& stream_out = std::cout) const {
    stream_out << "type: \t\t\t" << DensityEstimationTypeParser::toString(type_)
               << std::endl;
    stream_out << "decomposition: \t"
               << datadriven::MatrixDecompositionTypeParser::toString(
                      decomposition_)
               << std::endl;
    stream_out << "useOfflinePermutation: \t" << std::boolalpha
               << useOfflinePermutation << std::endl;
    stream_out << "normalize: \t\t" << std::boolalpha << normalize_
               << std::endl;
    stream_out << "iCholSweepsDecompose \t" << iCholSweepsDecompose_
               << std::endl;
    stream_out << "iCholSweepsRefine \t" << iCholSweepsRefine_ << std::endl;
    stream_out << "iCholSweepsUpdateLambda " << iCholSweepsUpdateLambda_
               << std::endl;
    stream_out << "iCholSweepsSolver \t" << iCholSweepsSolver_ << std::endl;
  }
  */
};

}  // namespace datadriven
}  // namespace sgpp

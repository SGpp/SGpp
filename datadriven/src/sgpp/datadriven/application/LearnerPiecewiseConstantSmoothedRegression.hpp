// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/application/DensityEstimator.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/datadriven/operation/hash/OperationPiecewiseConstantRegression/Node.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {

// --------------------------------------------------------------------------

class LearnerPiecewiseConstantSmoothedRegression {
 private:
  sgpp::base::RegularGridConfiguration gridConfig;

  sgpp::base::AdpativityConfiguration adaptivityConfig;

  sgpp::solver::SLESolverConfiguration solverConfig;

  sgpp::datadriven::RegularizationConfiguration regularizationConfig;

  bool verbose;

 public:
  /**
   * Constructor
   *
   * @param gridConfig grid configuration
   * @param adaptivityConfig adaptive refinement configuration
   * @param solverConfig solver configuration (CG)
   * @param regularizationConfig config for regularization operator
   * @param verbose report additional information on the console
   */
  LearnerPiecewiseConstantSmoothedRegression(sgpp::base::RegularGridConfiguration&
      gridConfig,
      sgpp::base::AdpativityConfiguration& adaptivityConfig,
      sgpp::solver::SLESolverConfiguration& solverConfig,
      sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
      bool verbose);

  /**
   * Does the learning step on a given grid, training set and regularization parameter lambda
   *
   * @param piecewiseRegressor
   * @param grid grid
   * @param alpha coefficient vector
   * @param lambda regularization parameter
   */
  void train(sgpp::datadriven::PiecewiseConstantRegression::Node&
             piecewiseRegressor, sgpp::base::Grid& grid,
             sgpp::base::DataVector& alpha, double lambda);

  /**
   * generates the regularization matrix
   * @param grid grid
   */
  sgpp::base::OperationMatrix* computeRegularizationMatrix(
    sgpp::base::Grid& grid);

  /**
   * Does cross-validation to obtain a suitable regularization parameter
   */
  double optimizeLambdaCV(size_t kFold);

  /**
   * splits the complete sample set in a set of smaller training and test
   * samples for cross-validation.
   *
   * @param strain vector containing the training samples for cv
   * @param stest vector containing the test samples for cv
   */
  void splitset(std::vector<sgpp::base::DataMatrix*>& strain,
                std::vector<sgpp::base::DataMatrix*>& stest);
};

}  // namespace datadriven
}  // namespace sgpp


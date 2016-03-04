// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LEARNERSGDELOG_HPP_
#define LEARNERSGDELOG_HPP_

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <sgpp/globaldef.hpp>

#include "DensityEstimator.hpp"
#include "LearnerSGDE.hpp"

namespace sgpp {
namespace datadriven {

class LearnerSGDElog : public datadriven::LearnerSGDE {
 public:
  /**
   * Constructor
   *
   * @param gridConfig grid configuration
   * @param adaptivityConfig adaptive refinement configuration
   * @param solverConfig solver configuration (CG)
   * @param regularizationConfig config for regularization operator
   * @param learnerSGDEConfig configuration for the learner
   */
  LearnerSGDElog(sgpp::base::RegularGridConfiguration& gridConfig,
                 sgpp::base::AdpativityConfiguration& adaptivityConfig,
                 sgpp::solver::SLESolverConfiguration& solverConfig,
                 sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
                 LearnerSGDEConfiguration& learnerSGDEConfig);

  virtual ~LearnerSGDElog();

  /**
   * This methods evaluates the sparse grid density at a single point
   * @param x DataVector length equal to dimensionality
   */
  virtual double pdf(sgpp::base::DataVector& x);

 protected:
  /**
   * Does the learning step on a given grid, training set and regularization parameter lambda
   *
   * @param grid grid
   * @param alpha coefficient vector
   * @param train sample set
   * @param lambdaReg regularization parameter
   */
  virtual void train(sgpp::base::Grid& grid, sgpp::base::DataVector& alpha,
                     sgpp::base::DataMatrix& train, double lambdaReg);

  double norm(sgpp::base::DataVector& v, sgpp::base::DataVector& w);
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* LEARNERSGDELOG_HPP_ */

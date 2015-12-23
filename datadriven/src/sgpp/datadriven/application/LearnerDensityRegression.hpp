// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include "DensityEstimator.hpp"
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/datadriven/operation/hash/OperationOcttreeHistogramRegression/Node.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/pde/application/RegularizationConfiguration.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

// --------------------------------------------------------------------------

class LearnerDensityRegression {
private:

    SGPP::base::RegularGridConfiguration gridConfig;

    SGPP::base::AdpativityConfiguration adaptivityConfig;

    SGPP::solver::SLESolverConfiguration solverConfig;

    SGPP::pde::RegularizationConfiguration regularizationConfig;

    bool verbose;
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
    LearnerDensityRegression(SGPP::base::RegularGridConfiguration& gridConfig,
    SGPP::base::AdpativityConfiguration& adaptivityConfig,
    SGPP::solver::SLESolverConfiguration& solverConfig,
    SGPP::pde::RegularizationConfiguration& regularizationConfig, bool verbose);

    /**
     * Does the learning step on a given grid, training set and regularization parameter lambda
     *
     * @param grid grid
     * @param alpha coefficient vector
     * @param train sample set
     * @param lambdaReg regularization parameter
     */
    void train(SGPP::datadriven::HistogramTree::Node &piecewiseRegressor, SGPP::base::Grid& grid,
    SGPP::base::DataVector& alpha, float_t lambda);

    /**
     * generates the regularization matrix
     * @param grid grid
     */
    SGPP::base::OperationMatrix* computeRegularizationMatrix(
    SGPP::base::Grid& grid);

};

} /* namespace datadriven */
} /* namespace SGPP */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/application/DensityEstimator.hpp>
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
     * @param verbose report additional information on the console
     */
    LearnerDensityRegression(SGPP::base::RegularGridConfiguration& gridConfig,
    SGPP::base::AdpativityConfiguration& adaptivityConfig,
    SGPP::solver::SLESolverConfiguration& solverConfig,
    SGPP::pde::RegularizationConfiguration& regularizationConfig, bool verbose);

    /**
     * Does the learning step on a given grid, training set and regularization parameter lambda
     *
     * @param piecewiseRegressor
     * @param grid grid
     * @param alpha coefficient vector
     * @param lambda regularization parameter
     */
    void train(SGPP::datadriven::HistogramTree::Node &piecewiseRegressor, SGPP::base::Grid& grid,
    SGPP::base::DataVector& alpha, float_t lambda);

    /**
     * generates the regularization matrix
     * @param grid grid
     */
    SGPP::base::OperationMatrix* computeRegularizationMatrix(
    SGPP::base::Grid& grid);

    /**
     * Does cross-validation to obtain a suitable regularization parameter
     */
    float_t optimizeLambdaCV(size_t kFold);

    /**
     * splits the complete sample set in a set of smaller training and test
     * samples for cross-validation.
     *
     * @param strain vector containing the training samples for cv
     * @param stest vector containing the test samples for cv
     */
    void splitset(std::vector<SGPP::base::DataMatrix*>& strain,
                  std::vector<SGPP::base::DataMatrix*>& stest);

};

} /* namespace datadriven */
} /* namespace SGPP */

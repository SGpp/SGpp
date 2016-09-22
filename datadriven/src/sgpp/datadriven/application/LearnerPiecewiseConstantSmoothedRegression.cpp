// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/LearnerPiecewiseConstantSmoothedRegression.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/hash/OperationFirstMoment.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/datadriven/algorithm/PiecewiseConstantSmoothedRegressionSystemMatrix.hpp>

#include <stddef.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

using sgpp::base::GridStorage;
using sgpp::base::GridGenerator;
using sgpp::base::DataVector;
using sgpp::base::OperationMatrix;
using sgpp::base::Grid;
using sgpp::base::SurplusRefinementFunctor;
using sgpp::base::GridPoint;
using sgpp::base::OperationEval;
using sgpp::base::application_exception;

using std::endl;
using std::cout;

namespace sgpp {
namespace datadriven {

LearnerPiecewiseConstantSmoothedRegression::LearnerPiecewiseConstantSmoothedRegression(
  sgpp::base::RegularGridConfiguration& gridConfig,
  sgpp::base::AdpativityConfiguration& adaptivityConfig,
  sgpp::solver::SLESolverConfiguration& solverConfig,
  sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
  bool verbose) :
  gridConfig(gridConfig), adaptivityConfig(adaptivityConfig),
  solverConfig(solverConfig), regularizationConfig(
    regularizationConfig), verbose(verbose) {
}

// ---------------------------------------------------------------------------

void LearnerPiecewiseConstantSmoothedRegression::train(
  sgpp::datadriven::PiecewiseConstantRegression::Node& piecewiseRegressor,
  Grid& grid,
  DataVector& alpha, double lambda) {
  size_t dim = grid.getDimension();

  GridStorage* gridStorage = &grid.getStorage();
  GridGenerator& gridGen = grid.getGenerator();
  DataVector rhs(grid.getSize());
  alpha.resize(grid.getSize());
  alpha.setAll(0.0);

  if (verbose) {
    cout << "# LearnerDensityRegression: grid points " << grid.getSize() << endl;
  }

  for (size_t ref = 0; ref <= adaptivityConfig.numRefinements_; ref++) {
    OperationMatrix* C = computeRegularizationMatrix(grid);

    sgpp::datadriven::PiecewiseConstantSmoothedRegressionSystemMatrix SMatrix(
      piecewiseRegressor, grid, *C, lambda);

    if (verbose) {
      cout << "# integrating rhs" << std::endl;
    }

    SMatrix.generateb(rhs);

    if (verbose) {
      cout << "# LearnerDensityRegression: Solving " << endl;
    }

    sgpp::solver::ConjugateGradients myCG(solverConfig.maxIterations_,
                                          solverConfig.eps_);
    myCG.solve(SMatrix, alpha, rhs, false, false, solverConfig.threshold_);

    if (ref < adaptivityConfig.numRefinements_) {
      if (verbose) {
        cout << "# LearnerDensityRegression: Refine grid ... " << std::endl;
      }

      // Weight surplus with function evaluation at grid points
      std::unique_ptr<OperationEval> opEval(sgpp::op_factory::createOperationEval(grid));
      DataVector p(dim);
      DataVector alphaWeight(alpha.getSize());

      for (size_t i = 0; i < gridStorage->getSize(); i++) {
        gridStorage->getPoint(i).getStandardCoordinates(p);
        alphaWeight[i] = alpha[i] * opEval->eval(alpha, p);
      }

      SurplusRefinementFunctor srf(alphaWeight, adaptivityConfig.noPoints_,
                                   adaptivityConfig.threshold_);
      gridGen.refine(srf);

      if (verbose) {
        cout << "# LearnerDensityRegression: ref " << ref << "/" <<
             adaptivityConfig.numRefinements_ - 1 << ": "
             << grid.getSize() << endl;
      }

      alpha.resize(grid.getSize());
      rhs.resize(grid.getSize());
      alpha.setAll(0.0);
      rhs.setAll(0.0);
    }

    delete C;
  }

  return;
}

OperationMatrix*
LearnerPiecewiseConstantSmoothedRegression::computeRegularizationMatrix(
  sgpp::base::Grid& grid) {
  OperationMatrix* C = NULL;

  if (regularizationConfig.regType_ ==
      sgpp::datadriven::RegularizationType::Identity) {
    C = sgpp::op_factory::createOperationIdentity(grid);
  } else if (regularizationConfig.regType_ ==
             sgpp::datadriven::RegularizationType::Laplace) {
    C = sgpp::op_factory::createOperationLaplace(grid);
  } else {
    throw application_exception("LearnerDensityRegression::train : unknown regularization type");
  }

  return C;
}

}  // namespace datadriven
}  // namespace sgpp

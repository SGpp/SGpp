// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

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
#include <sgpp/datadriven/algorithm/LogDensitySystemMatrix.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>

#include "LearnerSGDElog.hpp"
#include "LearnerSGDE.hpp"

namespace sgpp {
namespace datadriven {

LearnerSGDElog::LearnerSGDElog(base::RegularGridConfiguration& gridConfig,
                               base::AdpativityConfiguration& adaptivityConfig,
                               solver::SLESolverConfiguration& solverConfig,
                               datadriven::RegularizationConfiguration& regularizationConfig,
                               CrossvalidationForRegularizationConfiguration& crossvalidationConfig)
    : datadriven::LearnerSGDE(gridConfig, adaptivityConfig, solverConfig, regularizationConfig,
                              crossvalidationConfig) {}

LearnerSGDElog::~LearnerSGDElog() {}

// -----------------------------------------------------------------------

double LearnerSGDElog::norm(base::DataVector& v, base::DataVector& w) {
  double res = 0.0;

  for (size_t i = 0; i < v.getSize(); i++) {
    res += (v[i] - w[i]) * (v[i] - w[i]);
  }

  return res;
}

void LearnerSGDElog::train(base::Grid& grid, base::DataVector& alpha, base::DataMatrix& data,
                           double lambdaReg) {
  alpha.resize(grid.getSize());
  alpha.setAll(0.0);

  if (!crossvalidationConfig.silent_) {
    std::cout << "# LearnerSGDE: grid points " << grid.getSize() << std::endl;
  }

  // learn the usual density
  LearnerSGDE::train(grid, alpha, data, lambdaReg);

  // learn now the log density
  std::unique_ptr<base::OperationMatrix> C = computeRegularizationMatrix(grid);
  base::DataVector rhs(grid.getSize());

  // ---------------------------------------------------
  // learning
  datadriven::LogDensitySystemMatrix SMatrix(grid, alpha, data, *C, 1e-5);
  SMatrix.generateb(rhs);
  solver::ConjugateGradients myCG(solverConfig.maxIterations_, solverConfig.eps_);
  myCG.solve(SMatrix, alpha, rhs, false, false, solverConfig.threshold_);
  // ---------------------------------------------------
}

double LearnerSGDElog::pdf(base::DataVector& x) {
  std::unique_ptr<base::OperationEval> opEval = op_factory::createOperationEval(*grid);
  double ret = std::exp(opEval->eval(*alpha, x));
  return ret;
}

} /* namespace datadriven */
} /* namespace sgpp */

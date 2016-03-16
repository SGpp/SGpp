// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/application/LearnerSGD.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

namespace sgpp {
namespace datadriven {

LearnerSGD::LearnerSGD(sgpp::datadriven::RegularizationType& regularization,
                       const bool isRegression, const bool isVerbose)
    : Learner(regularization, isRegression, isVerbose) {}

void LearnerSGD::train(sgpp::base::DataMatrix& trainDataset, sgpp::base::DataVector& classes,
                       sgpp::base::RegularGridConfiguration& GridConfig, size_t maxIterations,
                       double eps, double lambda, double gamma) {
  // using namespace sgpp::base;
  using sgpp::base::DataVector;

  // Initialize Grid
  InitializeGrid(GridConfig);

  if (grid == NULL) return;

  alpha->setAll(0.0);

  size_t num_coeff = grid->getSize();
  size_t dim = trainDataset.getNcols();

  std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEval(*grid));

  // Execute SGD
  size_t numIterations = 0;

  while (numIterations < maxIterations) {
    // Get random x and y pair
    int k = getRandom(static_cast<int>(trainDataset.getNrows()) - 1);
    DataVector x(dim);
    trainDataset.getRow((size_t)k, x);
    double y = classes[k];

    // Calculate delta^n according to [Maier BA, 5.10]:

    // tmp1 = (b_k^T * a^n - y_k) where
    // b_k = (phi_1(x_k) ... phi_N(x_k))
    double tmp1 = opEval->eval(*alpha, x) - y;

    // delta^n = 2 * gamma * (b_k * tmp1 + lambda * a^n)
    DataVector delta(num_coeff);

    for (unsigned int i = 0; i < num_coeff; i++) {
      DataVector unit_alpha(num_coeff);
      unit_alpha.setAll(0.0);
      unit_alpha[i] = 1;

      delta[i] = 2 * gamma * (opEval->eval(unit_alpha, x) * tmp1 + lambda * (*alpha)[i]);
    }

    // update alpha
    // a^{n+1} = a^n - delta^n
    for (unsigned int i = 0; i < num_coeff; i++) {
      (*alpha)[i] = (*alpha)[i] - delta[i];
    }

    // check if below eps
    bool is_below_eps = true;

    for (unsigned int i = 0; i < num_coeff; i++) {
      if (fabs(delta[i]) >= eps) {
        is_below_eps = false;
        break;
      }
    }

    if (is_below_eps) {
      return;
    }

    numIterations++;
  }

  isTrained = true;
}

int LearnerSGD::getRandom(int limit) {
  int divisor = RAND_MAX / (limit + 1);
  int r;

  do {
    r = rand() / divisor;
  } while (r > limit);

  return r;
}

sgpp::base::DataVector& LearnerSGD::getAlpha() { return *alpha; }

sgpp::base::Grid& LearnerSGD::getGrid() { return *grid; }

LearnerSGD::~LearnerSGD() {}

}  // namespace datadriven
}  // namespace sgpp

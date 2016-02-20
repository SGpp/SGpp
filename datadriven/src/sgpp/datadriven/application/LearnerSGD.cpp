// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/application/LearnerSGD.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>


namespace SGPP {
namespace datadriven {

LearnerSGD::LearnerSGD(SGPP::datadriven::RegularizationType& regularization,
                       const bool isRegression, const bool isVerbose) : Learner(regularization,
                             isRegression, isVerbose) {
}

void LearnerSGD::train(
    SGPP::base::DataMatrix& trainDataset,
    SGPP::base::DataVector& classes,
    SGPP::base::RegularGridConfiguration& GridConfig,
    size_t maxIterations,
    float_t eps,
    float_t lambda,
    float_t gamma) {
  // using namespace SGPP::base;
  using SGPP::base::DataVector;

  // Initialize Grid
  InitializeGrid(GridConfig);

  if (grid_ == NULL)
    return;

  alpha_->setAll(0.0);

  size_t num_coeff = grid_->getStorage().size();
  size_t dim = trainDataset.getNcols();

  std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEval(*grid_));

  // Execute SGD
  size_t numIterations = 0;

  while (numIterations < maxIterations) {
    // Get random x and y pair
    int k = getRandom(static_cast<int>(trainDataset.getNrows()) - 1);
    DataVector x(dim);
    trainDataset.getRow((size_t)k, x);
    float_t y = classes[k];

    // Calculate delta^n according to [Maier BA, 5.10]:

    // tmp1 = (b_k^T * a^n - y_k) where
    // b_k = (phi_1(x_k) ... phi_N(x_k))
    float_t tmp1 = opEval->eval(*alpha_, x) - y;

    // delta^n = 2 * gamma * (b_k * tmp1 + lambda * a^n)
    DataVector delta(num_coeff);

    for (unsigned int i = 0; i < num_coeff; i++) {
      DataVector unit_alpha(num_coeff);
      unit_alpha.setAll(0.0);
      unit_alpha[i] = 1;

      delta[i] = 2 * gamma * (opEval->eval(unit_alpha, x) * tmp1 + lambda * (*alpha_)[i]);
    }

    // update alpha
    // a^{n+1} = a^n - delta^n
    for (unsigned int i = 0; i < num_coeff; i++) {
      (*alpha_)[i] = (*alpha_)[i] - delta[i];
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

  isTrained_ = true;
}

int LearnerSGD::getRandom(int limit) {
  int divisor = RAND_MAX / (limit + 1);
  int r;

  do {
    r = rand() / divisor;
  } while (r > limit);

  return r;
}

SGPP::base::DataVector* LearnerSGD::getAlpha() {
  return alpha_;
}

SGPP::base::Grid* LearnerSGD::getGrid() {
  return grid_;
}

LearnerSGD::~LearnerSGD() {
}

}  // namespace datadriven
}  // namespace SGPP


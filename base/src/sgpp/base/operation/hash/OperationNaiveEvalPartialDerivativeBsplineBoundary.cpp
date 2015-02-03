// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalPartialDerivativeBsplineBoundary.hpp>

namespace SGPP {
  namespace base {

    double OperationNaiveEvalPartialDerivativeBsplineBoundary::evalPartialDerivative(
      DataVector& alpha, const std::vector<double>& point, size_t deriv_dim) {
      const size_t n = storage->size();
      const size_t d = storage->dim();
      double result = 0.0;

      for (size_t i = 0; i < n; i++) {
        const GridIndex* gp = storage->get(i);
        double cur_val = 1.0;

        for (size_t t = 0; t < d; t++) {
          double val1d = ((t == deriv_dim) ?
                          base.evalDx(gp->getLevel(t), gp->getIndex(t), point[t]) :
                          base.eval(gp->getLevel(t), gp->getIndex(t), point[t]));

          if (val1d == 0.0) {
            cur_val = 0.0;
            break;
          }

          cur_val *= val1d;
        }

        result += alpha[i] * cur_val;
      }

      return result;
    }

  }
}

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalGradientBsplineClenshawCurtis.hpp>

namespace SGPP {
  namespace base {

    double OperationNaiveEvalGradientBsplineClenshawCurtis::evalGradient(
      DataVector& alpha, const std::vector<double>& point, DataVector& gradient) {
      const size_t n = storage->size();
      const size_t d = storage->dim();
      double result = 0.0;

      gradient.resize(storage->dim());
      gradient.setAll(0.0);

      DataVector cur_gradient(d);

      for (size_t i = 0; i < n; i++) {
        const GridIndex* gp = storage->get(i);
        double cur_val = 1.0;
        cur_gradient.setAll(alpha[i]);

        for (size_t t = 0; t < d; t++) {
          double val1d = base.eval(gp->getLevel(t), gp->getIndex(t), point[t]);
          double dx1d = base.evalDx(gp->getLevel(t), gp->getIndex(t), point[t]);

          cur_val *= val1d;

          for (size_t t2 = 0; t2 < d; t2++) {
            if (t2 == t) {
              cur_gradient[t2] *= dx1d;
            } else {
              cur_gradient[t2] *= val1d;
            }
          }
        }

        result += alpha[i] * cur_val;
        gradient.add(cur_gradient);
      }

      return result;
    }

  }
}

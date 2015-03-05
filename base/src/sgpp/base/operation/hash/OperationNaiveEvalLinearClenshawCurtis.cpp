// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalLinearClenshawCurtis.hpp>

namespace SGPP {
  namespace base {

    float_t OperationNaiveEvalLinearClenshawCurtis::eval(
      DataVector& alpha, std::vector<float_t>& point) {
      const size_t n = storage->size();
      const size_t d = storage->dim();
      float_t result = 0.0;

      for (size_t i = 0; i < n; i++) {
        const GridIndex* gp = storage->get(i);
        float_t curVal = 1.0;

        for (size_t t = 0; t < d; t++) {
          float_t val1d = base.eval(gp->getLevel(t), gp->getIndex(t), point[t]);

          if (val1d == 0.0) {
            curVal = 0.0;
            break;
          }

          curVal *= val1d;
        }

        result += alpha[i] * curVal;
      }

      return result;
    }

  }
}

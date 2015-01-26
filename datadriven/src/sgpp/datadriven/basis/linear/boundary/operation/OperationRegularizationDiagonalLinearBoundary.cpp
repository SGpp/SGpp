// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/datadriven/basis/linear/boundary/operation/OperationRegularizationDiagonalLinearBoundary.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <math.h>
//#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    OperationRegularizationDiagonalLinearBoundary::OperationRegularizationDiagonalLinearBoundary(base::GridStorage* storage, int mode, double k)
      : OperationRegularizationDiagonal(storage, mode, k) {
      init();
    }

    void OperationRegularizationDiagonalLinearBoundary::initHkmix (double k) {
      size_t dim = storage->dim();
      base::GridIndex* gi;
      double res;

      for (size_t i = 0; i < size; i++) {
        gi = storage->get(i);
        res = 1.0;

        for (size_t d = 0; d < dim; d++) {
          res *= pow(2, (2 * k - 1) * gi->getLevel(d) - 1);
        }

        diagonal[i] = res;
      }

    }

    void OperationRegularizationDiagonalLinearBoundary::initH0HkLaplace (double k) {
      size_t dim = storage->dim();
      base::GridIndex* gi;
      double res, resd;

      for (size_t i = 0; i < size; i++) {
        gi = storage->get(i);
        res = 0.0;

        for (size_t d = 0; d < dim; d++) {
          // Hk in dimension d
          resd = pow(2, (2 * k - 1) * gi->getLevel(d) - 1);

          // "H0" in remaining dimensions
          for (size_t d2 = 0; d2 < d; d2++) {
            resd *= pow(2.0, -1 - gi->getLevel(d2));
          }

          for (size_t d2 = d + 1; d2 < dim; d2++) {
            resd *= pow(2.0, -1 - gi->getLevel(d2));
          }

          res += resd;
        }

        diagonal[i] = res;
      }

    }

  }
}
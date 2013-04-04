/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

#include "datadriven/operation/OperationRegularizationDiagonal.hpp"
#include "base/exception/generation_exception.hpp"
#include <algorithm>

namespace sg {
  namespace datadriven {

    OperationRegularizationDiagonal::OperationRegularizationDiagonal(base::GridStorage* storage, int mode, double k) : mode(mode), k(k), size(storage->size()), storage(storage), diagonal(storage->size()) {

      // remember size of grid to check for changes in grid
      this->size = storage->size();
    }

    void OperationRegularizationDiagonal::mult(base::DataVector& alpha, base::DataVector& result) {
      if (size != storage->size()) {
        diagonal = base::DataVector(size);
        init();
      }

      // apply diagonal
      result.copyFrom(alpha);
      result.componentwise_mult(diagonal);

    }

    void OperationRegularizationDiagonal::init() {
      if (mode == HKMIX) {
        initHkmix(k);
      } else if (mode == H0HKLAPLACE) {
        initH0HkLaplace(k);
      } else if (mode == ISOTROPIC_PENALTY) {
        initIsotropicPenalty();
      } else if (mode == ANISOTROPIC_PENALTY) {
        initAnisotropicPenalty();
      } else {
        throw new base::generation_exception("OperationRegularizationDiagonal: Unknown mode specified!");
      }
    }

    void OperationRegularizationDiagonal::initIsotropicPenalty() {
      // \frac{1}{\max\{l_1,\dots,l_d\}-\min\{l_1,\dots,l_d\}+1}d
      size_t dim = storage->dim();
      base::GridIndex* gi;

      for (size_t i = 0; i < size; i++) {
        gi = storage->get(i);
        diagonal[i] = 1.0 / (gi->getLevelMax() - gi->getLevelMin() + 1) * (double)dim;
      }
    }
    void OperationRegularizationDiagonal::initAnisotropicPenalty() {
      // \frac{1}{2}\log(1+(\frac{\max\{l_1,\dots,l_d\}}{\max\{\min\{l_1,\dots,l_d\},1\}})d)
      size_t dim = storage->dim();
      base::GridIndex* gi;

      for (size_t i = 0; i < size; i++) {
        gi = storage->get(i);
        diagonal[i] = 0.5 * log(1. + gi->getLevelMax() / std::max(static_cast<int>(gi->getLevelMin()), 1) * (double)dim);
      }
    }

  }
}


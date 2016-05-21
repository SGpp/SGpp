// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationRegularizationDiagonal.hpp>
#include <sgpp/base/exception/generation_exception.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>

namespace sgpp {
namespace datadriven {

OperationRegularizationDiagonal::OperationRegularizationDiagonal(base::GridStorage* storage,
                                                                 int mode, double k)
    : mode(mode), k(k), size(storage->getSize()), storage(storage), diagonal(storage->getSize()) {
  // remember size of grid to check for changes in grid
  this->size = storage->getSize();
}

void OperationRegularizationDiagonal::mult(base::DataVector& alpha, base::DataVector& result) {
  if (size != storage->getSize()) {
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
    throw base::generation_exception(
        "OperationRegularizationDiagonal: Unknown mode specified!");
  }
}

void OperationRegularizationDiagonal::initIsotropicPenalty() {
  // \frac{1}{\max\{l_1,\dots,l_d\}-\min\{l_1,\dots,l_d\}+1}d
  size_t dim = storage->getDimension();

  for (size_t i = 0; i < size; i++) {
    base::GridIndex& gi = storage->getGridIndex(i);
    diagonal[i] = 1.0 / (gi.getLevelMax() - gi.getLevelMin() + 1) * static_cast<double>(dim);
  }
}
void OperationRegularizationDiagonal::initAnisotropicPenalty() {
  // \frac{1}{2}\log(1+(\frac{\max\{l_1,\dots,l_d\}}{\max\{\min\{l_1,\dots,l_d\},1\}})d)
  size_t dim = storage->getDimension();

  for (size_t i = 0; i < size; i++) {
    base::GridIndex& gi = storage->getGridIndex(i);
    diagonal[i] =
        0.5 * log(1. +
                  static_cast<double>(gi.getLevelMax()) /
                      static_cast<double>(std::max(static_cast<int>(gi.getLevelMin()), 1)) *
                      static_cast<double>(dim));
  }
}
}  // namespace datadriven
}  // namespace sgpp

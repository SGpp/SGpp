// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/combigrid/operation/OperationPoleHierarchisationLinear.hpp>

namespace sgpp {
namespace combigrid {

OperationPoleHierarchisationLinear::OperationPoleHierarchisationLinear() {
}

OperationPoleHierarchisationLinear::~OperationPoleHierarchisationLinear() {
}

void OperationPoleHierarchisationLinear::apply(base::DataVector& values, size_t start, size_t step,
    size_t count, level_t level, bool hasBoundary) {
  index_t hInv = static_cast<index_t>(1) << level;
  index_t h = 1;

  for (level_t l = level; l > 0; l--) {
    size_t k = start + step * (h - (hasBoundary ? 0 : 1));

    for (index_t i = 1; i < hInv; i += 2) {
      values[k] -= (values[k - step * h] + values[k + step * h]) / 2.0;
      k += 2 * step * h;
    }

    hInv /= 2;
    h *= 2;
  }
}

}  // namespace combigrid
}  // namespace sgpp

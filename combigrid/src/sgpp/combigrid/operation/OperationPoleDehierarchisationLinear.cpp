// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/combigrid/operation/OperationPoleDehierarchisationLinear.hpp>

namespace sgpp {
namespace combigrid {

OperationPoleDehierarchisationLinear::OperationPoleDehierarchisationLinear() {
}

OperationPoleDehierarchisationLinear::~OperationPoleDehierarchisationLinear() {
}

void OperationPoleDehierarchisationLinear::apply(base::DataVector& values, size_t start,
    size_t step, size_t count, level_t level, bool hasBoundary) {
  index_t hInv = 2;
  index_t h = static_cast<index_t>(1) << ((level > 0) ? (level - 1) : 0);

  for (level_t l = 1; l <= level; l++) {
    size_t k = start + step * (h - (hasBoundary ? 0 : 1));

    for (index_t i = 1; i < hInv; i += 2) {
      values[k] += (values[k - step * h] + values[k + step * h]) / 2.0;
      k += 2 * step * h;
    }

    hInv *= 2;
    h /= 2;
  }
}

}  // namespace combigrid
}  // namespace sgpp

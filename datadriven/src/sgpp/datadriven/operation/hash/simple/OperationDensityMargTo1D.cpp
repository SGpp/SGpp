// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationDensityMargTo1D.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {

void OperationDensityMargTo1D::computeMarginalizationIndices(std::vector<size_t>& dim_x,
                                                             size_t numDims,
                                                             std::vector<size_t>& margDims) {
  size_t count = 0;

  for (size_t idim = 0; idim < numDims; idim++) {
    size_t i = 0;
    bool marginalize = true;
    while (i < dim_x.size() && marginalize) {
      marginalize &= dim_x[i] != idim;
      i++;
    }

    if (marginalize) {
      margDims.push_back(idim - count);
      count++;
    }
  }
}
}  // namespace datadriven
}  // namespace sgpp

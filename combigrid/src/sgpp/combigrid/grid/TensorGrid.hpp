// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/definitions.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class TensorGrid {
  std::vector<base::DataVector> onedimGrids;

 public:
  explicit TensorGrid(std::vector<base::DataVector> const &onedimGrids)
      : onedimGrids(onedimGrids) {}

  base::DataVector getGridPoint(MultiIndex const &index) {
    base::DataVector result(onedimGrids.size());  // TODO(holzmudd): dimensionality check?

    for (size_t i = 0; i < onedimGrids.size(); ++i) {
      result[i] = onedimGrids[i][index[i]];
    }

    return result;
  }

  size_t numPointsInDim(size_t dim) { return onedimGrids[dim].getSize(); }

  MultiIndex numPoints() {
    MultiIndex result(onedimGrids.size());

    for (size_t i = 0; i < onedimGrids.size(); ++i) {
      result[i] = onedimGrids[i].getSize();
    }

    return result;
  }

  std::vector<base::DataVector> const &get1DGrids() { return onedimGrids; }
};

} /* namespace combigrid */
} /* namespace sgpp */

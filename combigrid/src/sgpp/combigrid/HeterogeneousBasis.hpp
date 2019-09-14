// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>

#include <memory>
#include <vector>

namespace sgpp {
namespace combigrid {

class HeterogeneousBasis {
 public:
  HeterogeneousBasis() : bases1d() {
  }

  explicit HeterogeneousBasis(
      const std::vector<std::unique_ptr<base::Basis<level_t, index_t>>>& bases1d) : bases1d() {
    for (const std::unique_ptr<base::Basis<level_t, index_t>>& basis1d : bases1d) {
      this->bases1d.push_back(basis1d.get());
    }
  }

  explicit HeterogeneousBasis(const std::vector<base::Basis<level_t, index_t>*>& bases1d) :
      bases1d(bases1d) {
  }

  HeterogeneousBasis(size_t dim, base::Basis<level_t, index_t>& basis1d) : bases1d(dim, &basis1d) {
  }

  inline double eval(LevelVector level, IndexVector index, base::DataVector point) const {
    double result = 1.0;
    const size_t dim = bases1d.size();

    for (size_t d = 0; d < dim; d++) {
      result *= bases1d[d]->eval(level[d], index[d], point[d]);
    }

    return result;
  }

 protected:
  std::vector<base::Basis<level_t, index_t>*> bases1d;
};

}  // namespace combigrid
}  // namespace sgpp

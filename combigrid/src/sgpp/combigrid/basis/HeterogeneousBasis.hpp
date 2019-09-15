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
      const std::vector<std::unique_ptr<base::Basis<level_t, index_t>>>& bases1d,
      bool isHierarchical = true) : bases1d(), isHierarchical_(isHierarchical) {
    for (const std::unique_ptr<base::Basis<level_t, index_t>>& basis1d : bases1d) {
      this->bases1d.push_back(basis1d.get());
    }
  }

  explicit HeterogeneousBasis(const std::vector<base::Basis<level_t, index_t>*>& bases1d,
      bool isHierarchical = true) : bases1d(bases1d), isHierarchical_(isHierarchical) {
  }

  HeterogeneousBasis(size_t dim, base::Basis<level_t, index_t>& basis1d,
      bool isHierarchical = true) : bases1d(dim, &basis1d), isHierarchical_(isHierarchical) {
  }

  inline static void hierarchizeLevelIndex(level_t& level, index_t& index) {
    if (index == 0) {
      level = 0;
    } else if (level > 0) {
      while (index % 2 == 0) {
        level--;
        index /= 2;
      }
    }
  }

  inline double eval(const LevelVector& level, const IndexVector& index,
      const base::DataVector& point) const {
    double result = 1.0;
    const size_t dim = bases1d.size();

    if (isHierarchical_) {
      for (size_t d = 0; d < dim; d++) {
        level_t l = level[d];
        index_t i = index[d];
        hierarchizeLevelIndex(l, i);
        result *= bases1d[d]->eval(l, i, point[d]);
      }
    } else {
      for (size_t d = 0; d < dim; d++) {
        result *= bases1d[d]->eval(level[d], index[d], point[d]);
      }
    }

    return result;
  }

  size_t getDimension() const {
    return bases1d.size();
  }

  const std::vector<base::Basis<level_t, index_t>*> getBases1d() const {
    return bases1d;
  }

  bool isHierarchical() const {
    return isHierarchical_;
  }

 protected:
  std::vector<base::Basis<level_t, index_t>*> bases1d;
  bool isHierarchical_;
};

}  // namespace combigrid
}  // namespace sgpp

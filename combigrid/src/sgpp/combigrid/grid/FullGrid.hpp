// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/basis/HeterogeneousBasis.hpp>

namespace sgpp {
namespace combigrid {

class FullGrid {
 public:
  FullGrid() : level(), hasBoundary_(true), basis() {
  }

  FullGrid(const LevelVector& level, const HeterogeneousBasis& basis, bool hasBoundary = true) :
      level(level), hasBoundary_(hasBoundary), basis(basis) {
  }

  const LevelVector& getLevel() const {
    return level;
  }

  LevelVector& getLevel() {
    return level;
  }

  void setLevel(const LevelVector& level) {
    this->level = level;
  }

  size_t getLevel(size_t d) const {
    return level[d];
  }

  index_t getMinIndex(size_t d) const {
    return (hasBoundary() ? 0 : 1);
  }

  void getMinIndex(IndexVector& index) const {
    for (size_t d = 0; d < level.size(); d++) {
      index[d] = getMinIndex(d);
    }
  }

  index_t getMaxIndex(size_t d) const {
    return (static_cast<index_t>(1) << level[d]) - getMinIndex(d);
  }

  void getMaxIndex(IndexVector& index) const {
    for (size_t d = 0; d < level.size(); d++) {
      index[d] = getMaxIndex(d);
    }
  }

  index_t getNumberOfIndexVectors(size_t d) const {
    return getMaxIndex(d) - getMinIndex(d) + 1;
  }

  index_t getNumberOfIndexVectors() const {
    index_t result = 1;

    for (size_t d = 0; d < level.size(); d++) {
      result *= getNumberOfIndexVectors(d);
    }

    return result;
  }

  size_t getDimension() const {
    return level.size();
  }

  bool hasBoundary() const {
    return hasBoundary_;
  }

  void setHasBoundary(bool hasBoundary) {
    this->hasBoundary_ = hasBoundary;
  }

  const HeterogeneousBasis& getBasis() const {
    return basis;
  }

  void setBasis(const HeterogeneousBasis& basis) {
    this->basis = basis;
  }

 protected:
  LevelVector level;
  bool hasBoundary_;
  HeterogeneousBasis basis;
};

}  // namespace combigrid
}  // namespace sgpp

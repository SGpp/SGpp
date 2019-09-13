// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/GeneralOperation.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>

namespace sgpp {
namespace combigrid {

class OperationPole : public GeneralOperation {
 public:
  OperationPole(level_t level, bool hasBoundary = true) :
      level(level), hasBoundary_(hasBoundary) {
  }

  ~OperationPole() override {
  }

  void apply(base::DataVector& values) override {
    apply(values, 0, 1, values.size());
  }

  void apply(base::DataVector& values, size_t start, size_t step, size_t count) override = 0;

  level_t getLevel() const {
    return level;
  }

  void setLevel(level_t level) {
    this->level = level;
  }

  bool hasBoundary() const {
    return hasBoundary_;
  }

  void setHasBoundary(bool hasBoundary) {
    this->hasBoundary_ = hasBoundary;
  }

 protected:
  level_t level;
  bool hasBoundary_;
};

}  // namespace combigrid
}  // namespace sgpp

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>

#include <memory>

namespace sgpp {
namespace combigrid {

class OperationPole {
 public:
  OperationPole() {
  }

  virtual ~OperationPole() {
  }

  virtual void apply(base::DataVector& values, level_t level, bool hasBoundary = true) {
    apply(values, 0, 1, values.size(), level, hasBoundary);
  }

  virtual void apply(base::DataVector& values, size_t start, size_t step, size_t count,
      level_t level, bool hasBoundary = true) = 0;
};

}  // namespace combigrid
}  // namespace sgpp

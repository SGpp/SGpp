// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/operation/OperationPole.hpp>

namespace sgpp {
namespace combigrid {

class OperationPoleHierarchisationLinear : public OperationPole {
 public:
  OperationPoleHierarchisationLinear();

  ~OperationPoleHierarchisationLinear() override;

  void apply(base::DataVector& values, size_t start, size_t step, size_t count,
      level_t level, bool hasBoundary = true) override;
};

}  // namespace combigrid
}  // namespace sgpp

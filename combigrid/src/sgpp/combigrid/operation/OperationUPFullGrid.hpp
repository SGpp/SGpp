// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/grid/FullGrid.hpp>
#include <sgpp/combigrid/operation/OperationPole.hpp>

#include <memory>
#include <vector>

namespace sgpp {
namespace combigrid {

class OperationUPFullGrid {
 public:
  OperationUPFullGrid(const FullGrid& grid,
      const std::vector<std::unique_ptr<OperationPole>>& operationPole);

  OperationUPFullGrid(const FullGrid& grid, const std::vector<OperationPole*>& operationPole);

  OperationUPFullGrid(const FullGrid& grid, OperationPole& operationPole);

  void apply(base::DataVector& values);

  const FullGrid& getGrid() const;

  void setGrid(const FullGrid& grid);

  const std::vector<OperationPole*>& getOperationPole() const;

  void setOperationPole(const std::vector<OperationPole*>& operationPole);

 protected:
  FullGrid grid;
  std::vector<OperationPole*> operationPole;
};

}  // namespace combigrid
}  // namespace sgpp

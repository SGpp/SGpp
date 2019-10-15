// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/grid/CombinationGrid.hpp>
#include <sgpp/combigrid/operation/OperationPole.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class OperationUPCombinationGrid {
 public:
  OperationUPCombinationGrid(const CombinationGrid& grid,
      const std::vector<std::unique_ptr<OperationPole>>& operationPole);

  OperationUPCombinationGrid(const CombinationGrid& grid,
      const std::vector<OperationPole*> operationPole);

  OperationUPCombinationGrid(const CombinationGrid& grid, OperationPole& operationPole);

  void apply(std::vector<base::DataVector>& values);

  const CombinationGrid& getGrid() const;

  void setGrid(const CombinationGrid& grid);

  const std::vector<OperationPole*>& getOperationPole() const;

  void setOperationPole(const std::vector<OperationPole*>& operationPole);

 protected:
  CombinationGrid grid;
  std::vector<OperationPole*> operationPole;
};

}  // namespace combigrid
}  // namespace sgpp

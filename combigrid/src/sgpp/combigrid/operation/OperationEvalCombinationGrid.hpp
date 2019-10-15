// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/grid/CombinationGrid.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class OperationEvalCombinationGrid {
 public:
  explicit OperationEvalCombinationGrid(const CombinationGrid& grid);

  double eval(const std::vector<base::DataVector>& surpluses,
      const base::DataVector& point);

  void eval(const std::vector<base::DataVector>& surpluses,
      const base::DataMatrix& points, base::DataVector& result);

  const CombinationGrid& getGrid() const;

  void setGrid(const CombinationGrid& grid);

 protected:
  CombinationGrid grid;
};

}  // namespace combigrid
}  // namespace sgpp

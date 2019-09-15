// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/combigrid/grid/FullGrid.hpp>

namespace sgpp {
namespace combigrid {

class OperationEvalFullGrid : public base::OperationEval {
 public:
  OperationEvalFullGrid();

  explicit OperationEvalFullGrid(const FullGrid& grid);

  ~OperationEvalFullGrid() override;

  double eval(const base::DataVector& surpluses, const base::DataVector& point) override;

  virtual void eval(const base::DataVector& surpluses, const base::DataMatrix& points,
      base::DataVector& result);

  const FullGrid& getGrid() const;

  void setGrid(const FullGrid& grid);

 protected:
  FullGrid grid;
};

}  // namespace combigrid
}  // namespace sgpp

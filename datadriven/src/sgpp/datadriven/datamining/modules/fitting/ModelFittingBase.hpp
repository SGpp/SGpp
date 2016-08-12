// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {

using sgpp::base::OperationMatrix;
using sgpp::base::Grid;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

namespace datadriven {

enum class FittingSolverState { refine, solve };

class ModelFittingBase {
 public:
  ModelFittingBase();

  virtual ~ModelFittingBase();

  /// new grid and new data set
  virtual void fit(Dataset& dataset) = 0;

  /// reuse the grid and assume old data set
  virtual void refine() = 0;

  /// reuse grid and new data set
  virtual void update(Dataset& dataset) = 0;

  /**
   *
   * @param sample
   * @return
   */
  virtual double evaluate(base::DataVector& sample);

  /**
   *
   * @param samples
   * @param result
   * @return
   */
  virtual void evaluate(base::DataMatrix& samples, base::DataVector& results) const;

  virtual const base::Grid& getGrid() const;
  virtual const base::DataVector& getSurpluses() const;

 protected:
  OperationMatrix* getRegularizationMatrix(RegularizationType regType);
  void initializeGrid(base::RegularGridConfiguration gridConfig);

  std::shared_ptr<base::Grid> grid;
  std::shared_ptr<base::DataVector> alpha;
};

} /* namespace datadriven */
} /* namespace sgpp */

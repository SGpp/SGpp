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
namespace datadriven {

class ModelFittingBase {
 public:
  ModelFittingBase();

  virtual ~ModelFittingBase();

  /// new grid and new data set
  virtual void fit(datadriven::Dataset& dataset) = 0;

  /// reuse the grid and assume old data set
  virtual void refine() = 0;

  /// reuse grid and new data set
  virtual void update(datadriven::Dataset& dataset) = 0;

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
  virtual void evaluate(sgpp::base::DataMatrix& samples, sgpp::base::DataVector& result);

  virtual std::shared_ptr<base::Grid> getGrid();
  virtual std::shared_ptr<base::DataVector> getSurpluses();

 protected:
  std::shared_ptr<base::OperationMatrix> getRegularizationMatrix(
      sgpp::datadriven::RegularizationType regType);
  void initializeGrid(base::RegularGridConfiguration gridConfig);

  std::shared_ptr<base::Grid> grid;
  std::shared_ptr<base::DataVector> alpha;
};

} /* namespace datadriven */
} /* namespace sgpp */

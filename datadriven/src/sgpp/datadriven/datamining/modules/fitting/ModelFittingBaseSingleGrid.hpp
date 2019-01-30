/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 * ModelFittingBaseSingleGrid.hpp
 *
 *  Created on: May 22, 2018
 *      Author: dominik
 */

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/solver/SLESolver.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <memory>

namespace sgpp {

using sgpp::base::Grid;
using sgpp::base::DataVector;

namespace datadriven {
/**
 * Base class for models operating on a single grid (i.e. Least-Squares-Regression and
 * density estimation
 */
class ModelFittingBaseSingleGrid : public ModelFittingBase {
 public:
  /**
   * Default constructor
   */
  ModelFittingBaseSingleGrid();

  // TODO(lettrich): fix this as soon as all member variables are copyable.
  /**
   * Copy constructor - we cannot deep copy all member variables yet.
   * @param rhs const reference to the scorer object to copy from.
   */
  ModelFittingBaseSingleGrid(const ModelFittingBaseSingleGrid& rhs) = delete;

  /**
   * Move constructor
   * @param rhs R-value reference to a scorer object to moved from.
   */
  ModelFittingBaseSingleGrid(ModelFittingBaseSingleGrid&& rhs) = default;

  // TODO(lettrich): fix this as soon as all member variables are copyable.
  /**
   * Copy assign operator - we cannot deep copy all member variables yet.
   * @param rhs const reference to the scorer object to copy from.
   * @return rerefernce to this with updated values.
   */
  ModelFittingBaseSingleGrid& operator=(const ModelFittingBaseSingleGrid& rhs) = delete;

  /**
   * Move assign operator
   * @param rhs R-value reference to an a scorer object to move from.
   * @return rerefernce to this with updated values.
   */
  ModelFittingBaseSingleGrid& operator=(ModelFittingBaseSingleGrid&& rhs) = default;

  /**
   * virtual destructor.
   */
  virtual ~ModelFittingBaseSingleGrid() = default;

  /**
   * Get the underlying grid object for the current model.
   * @return the grid object.
   */
  Grid& getGrid();

  /**
   * Get the surpluses of the current grid
   * @return vector of surpluses.
   */
  DataVector& getSurpluses();

 protected:
  /**
   * the sparse grid that approximates the data.
   */
  std::unique_ptr<Grid> grid;

  /**
   * hierarchical surpluses of the #grid.
   */
  DataVector alpha;
};

} /* namespace datadriven */
} /* namespace sgpp */

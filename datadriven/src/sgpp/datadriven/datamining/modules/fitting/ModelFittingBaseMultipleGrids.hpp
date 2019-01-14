// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBaseSingleGrid.hpp>

#include <vector>

using sgpp::base::Grid;
using sgpp::base::DataVector;

namespace sgpp {
namespace datadriven{
/**
 * Base Class for machine learning models working with multiple Grids, like used in the sparse grid
 * combination technique.
 */

class ModelFittingBaseMultipleGrids : public ModelFittingBase{
public:
	/**
	 * Default constructor
	 */
	ModelFittingBaseMultipleGrids();

	/**
	 * Copy constructor
	 */
	ModelFittingBaseMultipleGrids(const ModelFittingBaseMultipleGrids& rhs) = delete;

	/**
	 * Move constructor
	 * @param rhs R-value reference to a scorer object to moved from.
	 */
	ModelFittingBaseMultipleGrids(ModelFittingBaseMultipleGrids&& rhs) = default;

	 /**
	  * Get the underlying grid object for the current model.
	  * @return the grid object.
	  */
	  Grid& getGrid(size_t index);

	/**
	 * Get surpluses of a certain grid
	 * @param index from the desired grid
	 * @return returns vector of surpluses
	 */
	DataVector& getSurpluses(size_t index);



protected:
	/**
	 * vector of single grid models that approximate the data
	 */
	std::vector<std::unique_ptr<ModelFittingBaseSingleGrid>> grids;

	/**
	 * vector of weights of the single grids used to evaluate them together
	 */
	std::vector<double> weights;
};
} //namespace datadriven
} //namespace sgpp

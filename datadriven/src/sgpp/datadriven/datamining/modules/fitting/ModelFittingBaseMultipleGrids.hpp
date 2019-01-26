// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/*
 * ModelFittingBaseMultipleGrids.cpp
 *
 *  Created on: Jan 14, 2019
 *      Author: nico
 */

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBaseSingleGrid.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridPoint.hpp>

#include <vector>

using sgpp::base::Grid;
using sgpp::base::GeneralGridConfiguration;
using sgpp::base::DataVector;
using sgpp::base::HashGridPoint;

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
	 * @return vector of surpluses
	 */
	DataVector& getSurpluses(size_t index);

	/**
	 * Evaluate combined grids on a defined point
	 * @param HashGridPoint to evaluate
	 * @return the weighted sum of the HashGridPoints value in the grids
	 */
	double evaluate(DataVector& point);



	/**
	 * Add a model to the models vector
	 * @param GeneralGridConfig
	 * @param weight of the new added model
	 * @return index of the new model to find it in the models vector
	 */
	std::unique_ptr<ModelFittingBaseSingleGrid> createNewModel(GeneralGridConfiguration config);

protected:
	/**
	 * vector of single grid models that approximate the data
	 */
	std::vector<std::unique_ptr<ModelFittingBaseSingleGrid>> models;

	/**
	 * vector of weights of the single models used to evaluate them together
	 */
	std::vector<double> weights;

	/**
	 * Dimension of the Models Grids. All must have the same number of dimensions
	 */
	size_t dim;
};
} //namespace datadriven
} //namespace sgpp

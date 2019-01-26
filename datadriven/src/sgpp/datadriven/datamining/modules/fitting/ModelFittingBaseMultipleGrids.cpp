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

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBaseMultipleGrids.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBaseSingleGrid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>

#include <vector>

using sgpp::base::Grid;
using sgpp::base::DataVector;

namespace sgpp {
namespace datadriven{

ModelFittingBaseMultipleGrids::ModelFittingBaseMultipleGrids()
	: ModelFittingBase(), models{std::vector<std::unique_ptr<ModelFittingBaseSingleGrid>>{}}, weights{}, dim{0} {}

Grid& ModelFittingBaseMultipleGrids::getGrid(size_t index){
	return models[index]->getGrid();
}

DataVector& ModelFittingBaseMultipleGrids::getSurpluses(size_t index){
	return models[index]->getSurpluses();
}

double ModelFittingBaseMultipleGrids::evaluate(DataVector& point){
	//Data Vector must have same dimensions as the models of the grid
	if(point.size() != dim){
		throw sgpp::base::application_exception("ModelFittingBaseMultipleGrids::evaluate: Evaluation point dimensions don't fit grid dimensions");
	}
	double value = 0;
	for(size_t i = 0; i < models.size(); i++){
		value += models[i]->evaluate(point) * weights[i];
	}
	return value;
}

} //namespace sgpp
} //namespace datadriven


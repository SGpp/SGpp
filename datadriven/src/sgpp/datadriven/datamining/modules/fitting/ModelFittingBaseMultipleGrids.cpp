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

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBaseMultipleGrids.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBaseSingleGrid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <vector>

using sgpp::base::Grid;
using sgpp::base::DataVector;

namespace sgpp {
namespace datadriven{

ModelFittingBaseMultipleGrids::ModelFittingBaseMultipleGrids()
	: ModelFittingBase(), grids{std::vector<std::unique_ptr<ModelFittingBaseSingleGrid>>{}}, weights{} {}

Grid& ModelFittingBaseMultipleGrids::getGrid(size_t index){
	return grids[index]->getGrid();
}

DataVector& ModelFittingBaseMultipleGrids::getSurpluses(size_t index){
	return grids[index]->getSurpluses();
}

} //namespace sgpp
} //namespace datadriven


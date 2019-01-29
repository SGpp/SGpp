/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ModelFittingDensityEstimationMultiGrid.cpp
 *
 *  Created on: Jan 17, 2019
 *      Author: nico
 */

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationMultiGrid.hpp>
#include <vector>

namespace sgpp {
namespace datadriven {

ModelFittingDensityEstimationMultiGrid::ModelFittingDensityEstimationMultiGrid() {}

void ModelFittingDensityEstimationMultiGrid::fit(DataMatrix& dataset) {
  for (auto& model : models) {
    model->fit(dataset);
  }
}

void ModelFittingDensityEstimationMultiGrid::update(DataMatrix& samples) {
  for (auto& model : models) {
    model->update(samples);
  }
}

}  // namespace datadriven
}  // namespace sgpp

/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ModelFittingDensityEstimationCombiGrid.cpp
 *
 *  Created on: Jan 17, 2019
 *      Author: nico
 */

#include "ModelFittingDensityEstimationCombiGrid.hpp"

#include <sgpp/datadriven/datamining/configuration/CombiConfigurator.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp>
#include <vector>

using std::vector;
using std::unique_ptr;

namespace sgpp {
namespace datadriven {

ModelFittingDensityEstimationCombiGrid::ModelFittingDensityEstimationCombiGrid() {}

ModelFittingDensityEstimationCombiGrid::ModelFittingDensityEstimationCombiGrid(
    FitterConfigurationDensityEstimation& config) {
  auto combiconfig = vector<combiConfig>();
  CombiConfigurator combiconfigurator;
  combiconfigurator.getStandardCombi(combiconfig, config.getGridConfig().dim_,
                                     config.getGridConfig().level_);
  auto configtemp = config;
  models = vector<unique_ptr<ModelFittingDensityEstimationOnOff>>(combiconfig.size());
  for (size_t i = 1; i < combiconfig.size(); i++) {
    configtemp.getGridConfig().levelVector_ = combiconfig[i].levels;
    *models[i] = ModelFittingDensityEstimationOnOff(configtemp);
    weights[i] = combiconfig[i].coef;
  }
}

void ModelFittingDensityEstimationCombiGrid::fit(DataMatrix& dataset) {
  for (auto& model : models) {
    model->fit(dataset);
  }
}

void ModelFittingDensityEstimationCombiGrid::fit(Dataset& dataset) {
  for (auto& model : models) {
    model->fit(dataset);
  }
}

bool ModelFittingDensityEstimationCombiGrid::refine() { throw "not ready"; }

void ModelFittingDensityEstimationCombiGrid::reset() { throw "not ready"; }

bool ModelFittingDensityEstimationCombiGrid::refine(size_t newNoPoints,
                                                    std::list<size_t>* deletedGridPoints) {
  throw "not ready";
}

double ModelFittingDensityEstimationCombiGrid::evaluate(const DataVector& sample) {
  throw "not ready";
}

void ModelFittingDensityEstimationCombiGrid::evaluate(DataMatrix& samples, DataVector& results) {
  throw "not ready";
}

void ModelFittingDensityEstimationCombiGrid::update(DataMatrix& samples) {
  for (auto& model : models) {
    model->update(samples);
  }
}

void ModelFittingDensityEstimationCombiGrid::update(Dataset& dataset) { throw "not ready"; }

}  // namespace datadriven
}  // namespace sgpp

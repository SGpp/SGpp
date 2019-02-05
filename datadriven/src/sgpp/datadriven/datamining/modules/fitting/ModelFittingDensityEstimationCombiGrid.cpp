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

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/datadriven/datamining/configuration/CombiConfigurator.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCG.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp>

#include <iostream>
#include <vector>

using std::vector;
using std::unique_ptr;
using std::cout;
using sgpp::base::application_exception;

namespace sgpp {
namespace datadriven {

ModelFittingDensityEstimationCombiGrid::ModelFittingDensityEstimationCombiGrid() {}

ModelFittingDensityEstimationCombiGrid::ModelFittingDensityEstimationCombiGrid(
    FitterConfigurationDensityEstimation& config)
    : ModelFittingDensityEstimation{} {
  cout << "Creating ModelFittingDensityEstimationCombiGrid: \n";
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));
  vector<combiConfig> combiconfig;
  CombiConfigurator combiconfigurator;
  combiconfigurator.getStandardCombi(combiconfig, config.getGridConfig().dim_,
                                     config.getGridConfig().level_);

  models = vector<unique_ptr<ModelFittingDensityEstimation>>(combiconfig.size());
  weights = vector<double>(combiconfig.size());
  auto configtemp = FitterConfigurationDensityEstimation(config);
  for (size_t i = 0; i < combiconfig.size(); i++) {
    configtemp.getGridConfig().levelVector_ = combiconfig.at(i).levels;
    models.at(i) = createNewModel(configtemp);
    cout << "Model created\n";
    weights.at(i) = combiconfig.at(i).coef;
  }
}

std::unique_ptr<ModelFittingDensityEstimation>
ModelFittingDensityEstimationCombiGrid::createNewModel(
    sgpp::datadriven::FitterConfigurationDensityEstimation& densityEstimationConfig) {
  std::cout << "Creating New Model: \n";
  switch (DensityEstimationType::Decomposition) {  //!!!
    case DensityEstimationType::CG: {
      return std::make_unique<ModelFittingDensityEstimationCG>(densityEstimationConfig);
    }
    case DensityEstimationType::Decomposition: {
      return std::make_unique<ModelFittingDensityEstimationOnOff>(densityEstimationConfig);
    }
    default: { throw base::application_exception("Unknown density estimation type"); }
  }
}

void ModelFittingDensityEstimationCombiGrid::addNewModel(combiConfig config) {
  throw application_exception(
      "ModelFittingDensityEstimationCombiGrid::addNewModel(combiConfig config): not ready jet");
}

void ModelFittingDensityEstimationCombiGrid::fit(DataMatrix& newDataset) {
  for (auto& model : models) {
    model->fit(newDataset);
  }
}

void ModelFittingDensityEstimationCombiGrid::fit(Dataset& newDataset) {
  dataset = &newDataset;
  for (auto& model : models) {
    model->fit(newDataset);
  }
}

// bool ModelFittingDensityEstimationCombiGrid::refine() { throw "not ready"; }

void ModelFittingDensityEstimationCombiGrid::reset() {
  throw application_exception("ModelFittingDensityEstimationCombiGrid::reset(): not ready jet");
}

bool ModelFittingDensityEstimationCombiGrid::refine(size_t newNoPoints,
                                                    std::list<size_t>* deletedGridPoints) {
  throw application_exception(
      "ModelFittingDensityEstimationCombiGrid::refine(size_t newNoPoints, std::list<size_t>* "
      "deletedGridPoints): not ready jet");
}

double ModelFittingDensityEstimationCombiGrid::evaluate(const DataVector& sample) {
  throw application_exception(
      "ModelFittingDensityEstimationCombiGrid::evaluate(const DataVector& sample): not ready jet");
}

void ModelFittingDensityEstimationCombiGrid::evaluate(DataMatrix& samples, DataVector& results) {
  throw application_exception(
      "ModelFittingDensityEstimationCombiGrid::evaluate(DataMatrix& samples, DataVector& results): "
      "not ready jet");
}

void ModelFittingDensityEstimationCombiGrid::update(DataMatrix& newDataset) {
  for (auto& model : models) {
    model->update(newDataset);
  }
}

void ModelFittingDensityEstimationCombiGrid::update(Dataset& newDataset) {
  dataset = &newDataset;
  for (auto& model : models) {
      model->update(newDataset);
    }
}

bool ModelFittingDensityEstimationCombiGrid::isRefinable() {
  throw application_exception(
      "ModelFittingDensityEstimationCombiGrid::isRefinable(): not ready jet");
}

}  // namespace datadriven
}  // namespace sgpp

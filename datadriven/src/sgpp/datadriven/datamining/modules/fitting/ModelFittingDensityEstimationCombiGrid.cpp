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
  cout << "Creating ModelFittingDensityEstimationCombiGrid with Level = :"
       << config.getGridConfig().level_ << " \n";
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));
  components = vector<unique_ptr<ModelFittingDensityEstimation>>(0);
}

void ModelFittingDensityEstimationCombiGrid::fit(Dataset& newDataset) {
  dataset = &newDataset;
  fit(newDataset.getData());
}

void ModelFittingDensityEstimationCombiGrid::fit(DataMatrix& newDataset) {
  configurator = CombiConfigurator();
  configurator.initAdaptiveScheme(newDataset.getNcols(), config->getGridConfig().level_);
  configurator.getCombiScheme(componentConfigs);

  cout << "Printing the componentConfig:\n";
  for (auto v : componentConfigs) {
    cout << "[";
    for (auto b : v.levels) {
      cout << b << " ";
    }
    cout << "] " << v.coef << "\n";
  }

  components = vector<unique_ptr<ModelFittingDensityEstimation>>(componentConfigs.size());

  for (size_t i = 0; i < componentConfigs.size(); i++) {
    FitterConfigurationDensityEstimation newFitterConfig{};
    newFitterConfig.setupDefaults();
    newFitterConfig.getDensityEstimationConfig().type_ = config->getDensityEstimationConfig().type_;
    newFitterConfig.getRefinementConfig().numRefinements_ = 0;
    newFitterConfig.getGridConfig().levelVector_.clear();
    for (auto v : componentConfigs.at(i).levels) {
      cout << v << "\n";
      newFitterConfig.getGridConfig().levelVector_.push_back(v);
    }
    newFitterConfig.getGridConfig().generalType_ = config->getGridConfig().generalType_;

    components.at(i) = createNewModel(newFitterConfig);
    cout << "Model created\n";
  }
  cout << "Fitting the component grids to the Dataset: \n";
  for (auto& model : components) {
    model->fit(newDataset);
  }
  cout << "Done fitting the component grids. \n";
}

void ModelFittingDensityEstimationCombiGrid::update(Dataset& newDataset) {
  if (components.empty()) {
    fit(newDataset);
  }
  dataset = &newDataset;
  for (auto& model : components) {
    model->update(newDataset);
  }
}

void ModelFittingDensityEstimationCombiGrid::update(DataMatrix& newDataset) {
  if (components.empty()) {
    cout << "components.size() == 0 --> fit(newDataset)\n";
    fit(newDataset);
  }
  cout << "updating components..\n";
  for (auto& model : components) {
    model->update(newDataset);
  }
  cout << "done\n";
}

double ModelFittingDensityEstimationCombiGrid::evaluate(const DataVector& sample) {
  double result = 0;
  for (size_t i = 0; i < components.size(); i++) {
    result += components.at(i)->evaluate(sample) * componentConfigs.at(i).coef;
  }
  return result;
}

void ModelFittingDensityEstimationCombiGrid::evaluate(DataMatrix& samples, DataVector& results) {
  auto temp = sgpp::base::DataVector(results.size(), 0);
  results.setAll(0);
  for (size_t i = 0; i < components.size(); i++) {
    temp.setAll(0);
    components.at(i)->evaluate(samples, temp);
    temp.mult(componentConfigs.at(i).coef);
    results.add(temp);
  }
}

bool ModelFittingDensityEstimationCombiGrid::refine() {
  cout << "refine() refinementsPerfomed: " << refinementsPerformed;
  if (refinementsPerformed < config->getRefinementConfig().numRefinements_) {
    throw application_exception("ModelFittingDensityEstimationCombiGrid::refine(): not ready jet");
  }
  return false;
}

bool ModelFittingDensityEstimationCombiGrid::refine(size_t newNoPoints,
                                                    std::list<size_t>* deletedGridPoints) {
  throw application_exception(
      "ModelFittingDensityEstimationCombiGrid::refine(size_t newNoPoints, std::list<size_t>* "
      "deletedGridPoints): not ready jet\n");
}

void ModelFittingDensityEstimationCombiGrid::reset() {
  components.clear();
  refinementsPerformed = 0;
}

std::unique_ptr<ModelFittingDensityEstimation>
ModelFittingDensityEstimationCombiGrid::createNewModel(
    sgpp::datadriven::FitterConfigurationDensityEstimation& densityEstimationConfig) {
  std::cout << "Creating New ";
  switch (densityEstimationConfig.getDensityEstimationConfig().type_) {
    case DensityEstimationType::CG: {
      std::cout << "ModelFittingDensityEstimationCG\n";
      return std::make_unique<ModelFittingDensityEstimationCG>(densityEstimationConfig);
    }
    case DensityEstimationType::Decomposition: {
      std::cout << "ModelFittingDensityEstimationOnOff\n";
      return std::make_unique<ModelFittingDensityEstimationOnOff>(densityEstimationConfig);
    }
    default: { throw base::application_exception("Unknown density estimation type"); }
  }
}

void ModelFittingDensityEstimationCombiGrid::addNewModel(combiConfig config) {
  throw application_exception(
      "ModelFittingDensityEstimationCombiGrid::addNewModel(combiConfig config): not ready jet\n");
}

bool ModelFittingDensityEstimationCombiGrid::isRefinable() {
  throw application_exception(
      "ModelFittingDensityEstimationCombiGrid::isRefinable(): not ready jet\n");
}

}  // namespace datadriven
}  // namespace sgpp

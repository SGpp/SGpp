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
  this->config->getRefinementConfig().numRefinements_ = 0;
  models = vector<unique_ptr<ModelFittingDensityEstimation>>(0);
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

void ModelFittingDensityEstimationCombiGrid::fit(DataMatrix& newDataset) {
  cout << "voidi ModelFittingDensityEstimationCombiGrid::fit(DataMatrix& newDataset) ";
  vector<combiConfig> combiconfig;
  cout << "BREAK 1";
  CombiConfigurator combiconfigurator;
  cout << "BREAK 2";
  // FitterConfigurationDensityEstimation configtemp{};
  cout << "BREAK 3";
  // configtemp.setupDefaults();
  cout << "BREAK 4";
  // configtemp.getGridConfig() = config->getGridConfig();
  cout << "BREAK 5";
  // auto gridConfig = config->getGridConfig();

  combiconfigurator.getStandardCombi(combiconfig, newDataset.getNcols(),
                                     config->getGridConfig().level_);
  cout << "BREAK 6";

  models = vector<unique_ptr<ModelFittingDensityEstimation>>(combiconfig.size());
  cout << "BREAK 7";
  weights = vector<double>(combiconfig.size());
  cout << "BREAK 8";

  for (size_t i = 0; i < combiconfig.size(); i++) {
    FitterConfigurationDensityEstimation newFitterConfig{};
    newFitterConfig.setupDefaults();
    newFitterConfig.getDensityEstimationConfig().type_ = config->getDensityEstimationConfig().type_;
    newFitterConfig.getRefinementConfig().numRefinements_ = 0;
    newFitterConfig.getGridConfig().levelVector_ = combiconfig.at(i).levels;
    newFitterConfig.getGridConfig().generalType_ = config->getGridConfig().generalType_;

    models.at(i) = createNewModel(newFitterConfig);
    cout << "Model created\n";
    weights.at(i) = combiconfig.at(i).coef;
  }
  cout << "BREAK 9";
  for (auto& model : models) {
    model->fit(newDataset);
  }
  cout << "BREAK 10";
}

void ModelFittingDensityEstimationCombiGrid::fit(Dataset& newDataset) {
  cout << "void ModelFittingDensityEstimationCombiGrid::fit(Dataset& newDataset\n)";
  dataset = &newDataset;
  fit(newDataset.getData());
}

bool ModelFittingDensityEstimationCombiGrid::refine() {
  if (refinementsPerformed < config->getRefinementConfig().numRefinements_) {
    throw application_exception("ModelFittingDensityEstimationCombiGrid::refine(): not ready jet");
  }
  return false;
}

void ModelFittingDensityEstimationCombiGrid::reset() {
  throw application_exception("ModelFittingDensityEstimationCombiGrid::reset(): not ready jet\n");
}

bool ModelFittingDensityEstimationCombiGrid::refine(size_t newNoPoints,
                                                    std::list<size_t>* deletedGridPoints) {
  throw application_exception(
      "ModelFittingDensityEstimationCombiGrid::refine(size_t newNoPoints, std::list<size_t>* "
      "deletedGridPoints): not ready jet\n");
}

double ModelFittingDensityEstimationCombiGrid::evaluate(const DataVector& sample) {
  cout << "double ModelFittingDensityEstimationCombiGrid::evaluate(const DataVector& sample)\n";
  double result = 0;
  for (size_t i = 0; i < models.size(); i++) {
    result += models.at(i)->evaluate(sample) * weights.at(i);
  }
  return result;
}

void ModelFittingDensityEstimationCombiGrid::evaluate(DataMatrix& samples, DataVector& results) {
  cout << "void ModelFittingDensityEstimationCombiGrid::evaluate(DataMatrix& samples, DataVector& "
          "results)\n";
  auto temp = sgpp::base::DataVector(results.size(), 0);
  results.setAll(0);
  for (size_t i = 0; i < models.size(); i++) {
    temp.setAll(0);
    models.at(i)->evaluate(samples, temp);
    temp.mult(weights.at(i));
    results.add(temp);
  }
}

void ModelFittingDensityEstimationCombiGrid::update(DataMatrix& newDataset) {
  cout << "void ModelFittingDensityEstimationCombiGrid::update(DataMatrix& newDataset)\n";
  if (models.size() == 0) {
    fit(newDataset);
  }
  for (auto& model : models) {
    model->update(newDataset);
  }
}

void ModelFittingDensityEstimationCombiGrid::update(Dataset& newDataset) {
  cout << "void ModelFittingDensityEstimationCombiGrid::update(Dataset& newDataset)\n";
  if (models.size() == 0) {
    fit(newDataset);
  }
  dataset = &newDataset;
  for (auto& model : models) {
    model->update(newDataset);
  }
}

bool ModelFittingDensityEstimationCombiGrid::isRefinable() {
  throw application_exception(
      "ModelFittingDensityEstimationCombiGrid::isRefinable(): not ready jet\n");
}

}  // namespace datadriven
}  // namespace sgpp

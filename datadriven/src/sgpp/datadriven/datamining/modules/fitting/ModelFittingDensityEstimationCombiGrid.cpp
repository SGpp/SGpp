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
      newFitterConfig.getGridConfig().levelVector_.push_back(v);
    }
    newFitterConfig.getGridConfig().generalType_ = config->getGridConfig().generalType_;

    components.at(i) = createNewModel(newFitterConfig);
    cout << "Model created\n";
  }
  cout << "Fitting the component grids to the Dataset: \n";
  for (size_t i = 0; i < components.size(); i++) {
    components.at(i)->fit(newDataset);
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
    refinementsPerformed++;
    double max = 0;
    size_t ind = 0;
    /*
     * Finding the sub grid with the greatest error.
     * \TODO Add different kinds of error estimation
     */
    for (int i = 0; i < components.size(); i++) {
      double now =
          components.at(i)->getSurpluses().l2Norm() / components.at(i)->getSurpluses().getSize();
      cout << "Error: " << components.at(i)->getSurpluses().l2Norm() << " / "
           << components.at(i)->getSurpluses().getSize() << " = " << now << std::endl;
      if (now > max) {
        if (configurator.isRefinable(componentConfigs.at(i))) {
          max = now;
          ind = i;
        }
      }
    }
    /*
     * Refining the chosen block
     */
    configurator.refineBlock(componentConfigs.at(ind));
    vector<combiConfig> newConfigs;
    configurator.getCombiScheme(newConfigs);
    /*
     * Actualizing coefficients and finding newly-added and newly-remoced components
     */
    vector<bool> toAdd(newConfigs.size(), 1);
    vector<bool> toRemove(componentConfigs.size(), 1);

    for (size_t i = 0; i < newConfigs.size(); i++) {
      for (size_t k = 0; k < componentConfigs.size(); k++) {
        if (newConfigs.at(i).levels == componentConfigs.at(k).levels) {
          toAdd.at(i) = 0;
          if (toRemove.size() - k - 1 >= 0) {
            toRemove.at(toRemove.size() - k - 1) = 0;
          }
          componentConfigs.at(k).coef = newConfigs.at(i).coef;
        }
      }
    }
    /*
     * Removing components
     */
    for (size_t i = 0; i < toRemove.size(); i++) {
      removeModel(toRemove.at(i));
    }

    /*
     * Adding new components
     */
    for (size_t i = 0; i < toAdd.size(); i++) {
      if (toAdd.at(i)) {
        addNewModel(newConfigs.at(i));
      }
    }
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

void ModelFittingDensityEstimationCombiGrid::addNewModel(combiConfig combiconfig) {
  FitterConfigurationDensityEstimation newFitterConfig{};
  newFitterConfig.setupDefaults();
  newFitterConfig.getDensityEstimationConfig().type_ =
      this->config->getDensityEstimationConfig().type_;
  newFitterConfig.getRefinementConfig().numRefinements_ = 0;
  newFitterConfig.getGridConfig().levelVector_.clear();
  for (auto v : combiconfig.levels) {
    newFitterConfig.getGridConfig().levelVector_.push_back(v);
  }
  newFitterConfig.getGridConfig().generalType_ = this->config->getGridConfig().generalType_;

  components.push_back(createNewModel(newFitterConfig));
  componentConfigs.push_back(combiconfig);
  components.back()->fit(dataset->getData());
}

void ModelFittingDensityEstimationCombiGrid::removeModel(const size_t indexRev) {
  componentConfigs.erase(componentConfigs.end() - indexRev);
  components.erase(components.end() - indexRev);
}

bool ModelFittingDensityEstimationCombiGrid::isRefinable() {
  throw application_exception(
      "ModelFittingDensityEstimationCombiGrid::isRefinable(): not ready jet\n");
}

}  // namespace datadriven
}  // namespace sgpp

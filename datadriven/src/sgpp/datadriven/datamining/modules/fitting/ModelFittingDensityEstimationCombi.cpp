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

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/datadriven/datamining/configuration/CombiConfigurator.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCG.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCombi.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp>

#include <iostream>
#include <list>
#include <vector>

using std::vector;
using std::unique_ptr;
using sgpp::base::application_exception;

namespace sgpp {
namespace datadriven {

ModelFittingDensityEstimationCombi::ModelFittingDensityEstimationCombi() {}

ModelFittingDensityEstimationCombi::ModelFittingDensityEstimationCombi(
    FitterConfigurationDensityEstimation& config)
    : ModelFittingDensityEstimation{} {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));
  components = vector<unique_ptr<ModelFittingDensityEstimation>>(0);
  fitted = vector<bool>(0);
}

void ModelFittingDensityEstimationCombi::fit(Dataset& newDataset) {
  dataset = &newDataset;
  fit(newDataset.getData());
}

void ModelFittingDensityEstimationCombi::fit(DataMatrix& newDataset) {
  configurator = CombiConfigurator();
  configurator.initAdaptiveScheme(newDataset.getNcols(), config->getGridConfig().level_);
  configurator.getCombiScheme(componentConfigs);
  components = vector<unique_ptr<ModelFittingDensityEstimation>>(componentConfigs.size());
  fitted = vector<bool>(componentConfigs.size());

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
    fitted.at(i) = 0;
  }
  for (size_t i = 0; i < components.size(); i++) {
    components.at(i)->fit(newDataset);
    fitted.at(i) = 1;
  }
}

void ModelFittingDensityEstimationCombi::update(Dataset& newDataset) {
  if (components.empty()) {
    fit(newDataset);
  } else {
    for (size_t i = 0; i < components.size(); i++) {
      if (fitted.at(i) != 1) {
        components.at(i)->fit(newDataset.getData());
        fitted.at(i) = 1;
      }
    }
  }
  size_t gridpoints = 0;
  for (size_t i = 0; i < components.size(); i++) {
    gridpoints += components.at(i)->getGrid().getSize();
  }
  std::cout << "Refinement: " << refinementsPerformed << " Sum of Gridpoints: " << gridpoints
            << std::endl;
  dataset = &newDataset;
}

void ModelFittingDensityEstimationCombi::update(DataMatrix& newDataset) {
  if (components.empty()) {
    fit(newDataset);
  } else {
    for (size_t i = 0; i < components.size(); i++) {
      if (fitted.at(i) != 1) {
        components.at(i)->fit(newDataset);
        fitted.at(i) = 1;
      }
    }
  }
  size_t gridpoints = 0;
  for (size_t i = 0; i < components.size(); i++) {
    if (fitted.at(i)) {
      gridpoints += components.at(i)->getGrid().getSize();
    }
  }
  std::cout << "Refinement: " << refinementsPerformed << "; Sum of Gridpoints: " << gridpoints
            << std::endl;
}

double ModelFittingDensityEstimationCombi::evaluate(const DataVector& sample) {
  double result = 0;
  for (size_t i = 0; i < components.size(); i++) {
    if (fitted.at(i)) {
      result +=
          components.at(i)->evaluate(sample) * static_cast<double>(componentConfigs.at(i).coef);
    }
  }
  return result;
}

void ModelFittingDensityEstimationCombi::evaluate(DataMatrix& samples, DataVector& results) {
  auto temp = sgpp::base::DataVector(results.size(), 0);
  results.setAll(0);
  for (size_t i = 0; i < components.size(); i++) {
    if (fitted.at(i)) {
      temp.setAll(0);
      components.at(i)->evaluate(samples, temp);
      temp.mult(static_cast<double>(componentConfigs.at(i).coef));
      results.add(temp);
    }
  }
}

bool ModelFittingDensityEstimationCombi::refine() {
  if (componentConfigs.size() != components.size()) {
    throw base::application_exception("componentsConfig.size() != components.size()");
  }
  if (components.size() == 0) {
    throw base::application_exception("components.size() == 0");
  }

  if (refinementsPerformed < config->getRefinementConfig().numRefinements_) {
    refinementsPerformed++;

    /*
     * Finding the sub grid with the greatest error.
     * \TODO Add different kinds of error estimation
     */
    double max = 0;
    size_t ind = 0;
    for (size_t i = 0; i < components.size(); i++) {
      double now = components.at(i)->getSurpluses().l2Norm() /
                   static_cast<double>(components.at(i)->getSurpluses().getSize());
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
    configurator.refineComponent(componentConfigs.at(ind));
    vector<combiConfig> newConfigs;
    configurator.getCombiScheme(newConfigs);
    /*
     * Actualizing coefficients and finding newly-added and newly-removed components
     */
    vector<bool> toAdd(newConfigs.size(), 1);
    vector<bool> toRemove(componentConfigs.size(), true);

    for (size_t i = 0; i < newConfigs.size(); i++) {
      for (size_t k = 0; k < componentConfigs.size(); k++) {
        if (newConfigs.at(i).levels == componentConfigs.at(k).levels) {
          toAdd.at(i) = 0;
          toRemove.at(k) = false;
          componentConfigs.at(k).coef = newConfigs.at(i).coef;
        }
      }
    }
    /*
     * Removing components
     */
    for (size_t i = toRemove.size(); i > 0; i--) {
      if (toRemove.at(i - 1)) {
        removeModel(i - 1);
      }
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

bool ModelFittingDensityEstimationCombi::refine(size_t newNoPoints,
                                                std::list<size_t>* deletedGridPoints) {
  throw application_exception(
      "ModelFittingDensityEstimationCombiGrid::refine(size_t newNoPoints, std::list<size_t>* "
      "deletedGridPoints): not ready jet\n");
}

void ModelFittingDensityEstimationCombi::reset() {
  components.clear();
  refinementsPerformed = 0;
}

std::unique_ptr<ModelFittingDensityEstimation> ModelFittingDensityEstimationCombi::createNewModel(
    sgpp::datadriven::FitterConfigurationDensityEstimation& densityEstimationConfig) {
  switch (densityEstimationConfig.getDensityEstimationConfig().type_) {
    case DensityEstimationType::CG: {
      return std::make_unique<ModelFittingDensityEstimationCG>(densityEstimationConfig);
    }
    case DensityEstimationType::Decomposition: {
      return std::make_unique<ModelFittingDensityEstimationOnOff>(densityEstimationConfig);
    }
    default: { throw base::application_exception("Unknown density estimation type"); }
  }
}

void ModelFittingDensityEstimationCombi::addNewModel(combiConfig combiconfig) {
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
  fitted.push_back(0);
  componentConfigs.push_back(combiconfig);
}

void ModelFittingDensityEstimationCombi::removeModel(const size_t ind) {
  componentConfigs.erase(componentConfigs.begin() + ind);
  components.erase(components.begin() + ind);
  fitted.erase(fitted.begin() + ind);
}

bool ModelFittingDensityEstimationCombi::isRefinable() {
  throw application_exception(
      "ModelFittingDensityEstimationCombiGrid::isRefinable(): not ready jet\n");
}

}  // namespace datadriven
}  // namespace sgpp

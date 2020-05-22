// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCombi.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/datadriven/algorithm/CombiScheme.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCG.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp>

#include <iostream>
#include <list>
#include <utility>
#include <vector>

using sgpp::base::application_exception;
using std::unique_ptr;
using std::vector;

namespace sgpp {
namespace datadriven {

ModelFittingDensityEstimationCombi::ModelFittingDensityEstimationCombi() {}

ModelFittingDensityEstimationCombi::ModelFittingDensityEstimationCombi(
    const FitterConfigurationDensityEstimation& config)
    : ModelFittingDensityEstimation{} {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));
  components = vector<unique_ptr<ModelFittingDensityEstimation>>(0);
  fitted = vector<bool>(0);
  // If no object store is passed but the offline permutation is configured and the decomposition
  // type allows offline permutation, an object store is instantiated
  if (config.getDensityEstimationConfig().useOfflinePermutation_ &&
      DBMatOfflinePermutable::PermutableDecompositions.find(
          config.getDensityEstimationConfig().decomposition_) !=
          DBMatOfflinePermutable::PermutableDecompositions.end()) {
    this->objectStore = std::make_shared<DBMatObjectStore>();
    this->hasObjectStore = true;
  } else {
    this->hasObjectStore = false;
  }
}

ModelFittingDensityEstimationCombi::ModelFittingDensityEstimationCombi(
    const FitterConfigurationDensityEstimation& config,
    std::shared_ptr<DBMatObjectStore> objectStore)
    : ModelFittingDensityEstimationCombi(config) {
  this->hasObjectStore = true;
  this->objectStore = objectStore;
}

void ModelFittingDensityEstimationCombi::fit(Dataset& newDataset) {
  dataset = &newDataset;
  fit(newDataset.getData());
}

void ModelFittingDensityEstimationCombi::fit(DataMatrix& newDataset) {
  scheme.initialize(newDataset.getNcols(), config->getGridConfig().level_);
  componentConfigs = scheme.getCombiScheme();
  components = vector<unique_ptr<ModelFittingDensityEstimation>>(componentConfigs.size());
  fitted = vector<bool>(componentConfigs.size());

  for (size_t i = 0; i < componentConfigs.size(); i++) {
    FitterConfigurationDensityEstimation newFitterConfig{};
    newFitterConfig.setupDefaults();
    newFitterConfig.getRegularizationConfig().lambda_ =
        this->getFitterConfiguration().getRegularizationConfig().lambda_;
    newFitterConfig.getGridConfig().generalType_ = sgpp::base::GeneralGridType::ComponentGrid;
    newFitterConfig.getDensityEstimationConfig().decomposition_ =
        this->getFitterConfiguration().getDensityEstimationConfig().decomposition_;
    newFitterConfig.getDensityEstimationConfig().type_ = config->getDensityEstimationConfig().type_;
    newFitterConfig.getRefinementConfig().numRefinements_ = 0;
    newFitterConfig.getGridConfig().levelVector_.clear();
    for (auto v : componentConfigs.at(i).first) {
      newFitterConfig.getGridConfig().levelVector_.push_back(v);
    }

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
          components.at(i)->evaluate(sample) * static_cast<double>(componentConfigs.at(i).second);
    }
  }
  return result;
}

void ModelFittingDensityEstimationCombi::evaluate(DataMatrix& samples, DataVector& results) {
  auto temp = sgpp::base::DataVector(results.size(), 0);
  results.setAll(0.);
  for (size_t i = 0; i < components.size(); i++) {
    if (fitted.at(i)) {
      temp.setAll(0.);
      components.at(i)->evaluate(samples, temp);
      temp.mult(static_cast<double>(componentConfigs.at(i).second));
      results.add(temp);
    }
  }
}

bool ModelFittingDensityEstimationCombi::adapt() {
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
     * TODO: Add different kinds of error estimation
     */
    double max = 0;
    size_t ind = 0;
    for (size_t i = 0; i < components.size(); i++) {
      double now = components.at(i)->getSurpluses().l2Norm() /
                   static_cast<double>(components.at(i)->getSurpluses().getSize());
      if (now > max) {
        if (scheme.isRefinable(componentConfigs.at(i).first)) {
          max = now;
          ind = i;
        }
      }
    }

    /*
     * Refining the chosen block
     */
    scheme.refineComponent(componentConfigs.at(ind).first);
    std::vector<std::pair<std::vector<size_t>, int>> newConfigs = scheme.getCombiScheme();

    /*
     * Actualizing coefficients and finding newly-added and newly-removed components
     */
    vector<bool> toAdd(newConfigs.size(), 1);
    vector<bool> toRemove(componentConfigs.size(), true);

    for (size_t i = 0; i < newConfigs.size(); i++) {
      for (size_t k = 0; k < componentConfigs.size(); k++) {
        if (newConfigs.at(i).first == componentConfigs.at(k).first) {
          toAdd.at(i) = 0;
          toRemove.at(k) = false;
          componentConfigs.at(k).second = newConfigs.at(i).second;
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

bool ModelFittingDensityEstimationCombi::adapt(size_t newNoPoints,
                                               std::vector<size_t>& deletedGridPoints) {
  throw application_exception(
      "ModelFittingDensityEstimationCombiGrid::refine(size_t newNoPoints, std::vector<size_t>& "
      "deletedGridPoints): not implemented yet\n");
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
      if (this->hasObjectStore) {
        return std::make_unique<ModelFittingDensityEstimationOnOff>(densityEstimationConfig,
                                                                    objectStore);
      } else {
        return std::make_unique<ModelFittingDensityEstimationOnOff>(densityEstimationConfig);
      }
    }
  }

  throw base::application_exception("Unknown density estimation type");
}

void ModelFittingDensityEstimationCombi::addNewModel(
    std::pair<std::vector<size_t>, int> combiconfig) {
  FitterConfigurationDensityEstimation newFitterConfig{};
  newFitterConfig.setupDefaults();
  newFitterConfig.getDensityEstimationConfig().type_ =
      this->config->getDensityEstimationConfig().type_;
  newFitterConfig.getRefinementConfig().numRefinements_ = 0;
  newFitterConfig.getGridConfig().levelVector_.clear();
  for (auto v : combiconfig.first) {
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
      "ModelFittingDensityEstimationCombiGrid::isRefinable(): not implemented yet\n");
}

}  // namespace datadriven
}  // namespace sgpp

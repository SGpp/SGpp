/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ModelFittingDensityEstimation.cpp
 *
 * Created on: Jan 02, 2018
 *     Author: Kilian Röhner
 */

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>

using sgpp::base::Grid;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::SurplusRefinementFunctor;

using sgpp::base::application_exception;

namespace sgpp {
namespace datadriven {

ModelFittingDensityEstimation::ModelFittingDensityEstimation(
    const FitterConfigurationDensityEstimation& config)
    : ModelFittingBase{}, refinementsPerformed{0}, offline{nullptr} {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));
  // TODO(roehner): enable loading of decomposition from files
  offline = std::unique_ptr<DBMatOffline>{DBMatOfflineFactory::buildOfflineObject(
      this->config->getGridConfig(),
      this->config->getRefinementConfig(),
      this->config->getRegularizationConfig(),
      this->config->getDensityEstimationConfig())};
  offline->buildMatrix();
  // grid = offline->getGridPtr();
  offline->decomposeMatrix();
  online = std::unique_ptr<DBMatOnlineDE>{DBMatOnlineDEFactory::buildDBMatOnlineDE(*offline)};
  alpha = online->getAlpha();
}

// TODO(lettrich): exceptions have to be thrown if not valid.
double ModelFittingDensityEstimation::evaluate(const DataVector& sample) const {
  return online->eval(sample);;
}

// TODO(lettrich): exceptions have to be thrown if not valid.
void ModelFittingDensityEstimation::evaluate(DataMatrix& samples, DataVector& results) {
    online->eval(samples, results);
}

void ModelFittingDensityEstimation::fit(Dataset& newDataset) {
  online->computeDensityFunction(newDataset.getData());
}

bool ModelFittingDensityEstimation::refine() {
  /*if (grid != nullptr && offline->isRefineable()) {
    if (refinementsPerformed < config->getRefinementConfig().numRefinements_) {
      // create refinement functor
      SurplusRefinementFunctor refinementFunctor(alpha, config->getRefinementConfig().noPoints_,
                                                 config->getRefinementConfig().threshold_);
      // refine grid
      auto oldNoPoints = grid->getSize();
      grid->getGenerator().refine(refinementFunctor);
      auto newNoPoints = grid->getSize();
      if (newNoPoints > oldNoPoints) {
        // Tell the SLE manager that the grid changed (for interal data structures)
        alpha.resizeZero(newNoPoints);
        
        std::list<size_t> deletedGridPoints;
        // TODO(roehner) enable coarsening
        online->updateSystemMatrixDecomposition(newNoPoints - oldNoPoints,
                                                deletedGridPoints, online->getBestLambda());
        refinementsPerformed++;
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }

  } else {
    throw application_exception(
        "ModelFittingDensityEstimation: Can't refine before initial grid is created");
    return false;
  }*/
  return false;
}

void ModelFittingDensityEstimation::update(Dataset& newDataset) {
  fit(newDataset);
}

void ModelFittingDensityEstimation::resetState() { refinementsPerformed = 0; }

}  // namespace datadriven
}  // namespace sgpp

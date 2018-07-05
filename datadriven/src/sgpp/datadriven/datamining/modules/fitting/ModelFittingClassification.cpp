/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ModelFittingClassification.cpp
 *
 *  Created on: Jul 1, 2018
 *      Author: dominik
 */

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/datadriven/functors/MultiSurplusRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/GridPointBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/MultipleClassRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClassification.hpp>
#include <sgpp/datadriven/datamining/configuration/RefinementFunctorTypeParser.hpp>

#include <string>
#include <vector>
#include <list>
#include <map>

using sgpp::base::Grid;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

using sgpp::base::application_exception;

using sgpp::base::RefinementFunctorType;

namespace sgpp {
namespace datadriven {

ModelFittingClassification::ModelFittingClassification(
    const FitterConfigurationClassification& config)
    : refinementsPerformed{0} {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));
}

double ModelFittingClassification::evaluate(const DataVector& sample) {
  if (models.size() == 0) {
    std::string errorMessage = "Prediction impossible! No models were trained!";
    throw application_exception(errorMessage.c_str());
  } else {
    double prediction = 0.0, maxDensity = 0.0;
    bool evaluatedModel = false;
    for (auto& p : classIdx) {
      double label = p.first;
      size_t idx = p.second;
      if (classNumberInstances[idx] == 0) {
         // The model for this class was not trained -> no prediction possible for this model
         continue;
      }
      double classConditionalDensity = models[idx]->evaluate(sample);
      if (!evaluatedModel || classConditionalDensity > maxDensity) {
        maxDensity = classConditionalDensity;
        prediction = label;
      }
        evaluatedModel = true;
      }
      return prediction;
  }
}

void ModelFittingClassification::evaluate(DataMatrix& samples, DataVector& results) {
#pragma omp parallel for
  for (size_t i = 0; i < samples.getNrows(); i++) {
    DataVector tmp(samples.getNcols());
    samples.getRow(i, tmp);
    results.set(i, evaluate(tmp));
  }
}

void ModelFittingClassification::fit(Dataset& newDataset) {
  dataset = &newDataset;

  // Split the dataset into classes
  DataVector tmp(newDataset.getDimension());
  std::map<double, DataMatrix*> classSamples;
  for (size_t i = 0; i < newDataset.getNumberInstances(); i++) {
    double label = newDataset.getTargets().get(i);
    if (classSamples.find(label) == classSamples.end()) {
      classSamples[label] = new DataMatrix(0, newDataset.getDimension());
    }
    newDataset.getData().getRow(i, tmp);
    classSamples.at(label)->appendRow(tmp);
  }

  // Update the models
  for (auto& p : classSamples) {
    size_t idx = labelToIdx(p.first);
    DataMatrix *samples = p.second;
    models[idx]->update(*samples);
    classNumberInstances[idx] += samples->getNrows();
    delete samples;
  }
}

size_t ModelFittingClassification::labelToIdx(double label) {
  if (classIdx.find(label) == classIdx.end()) {
    // New class
    size_t idx = classIdx.size();
    classIdx[label] = idx;

    // Create a new model
    std::unique_ptr<ModelFittingDensityEstimation> model =
        std::make_unique<ModelFittingDensityEstimation>(
            dynamic_cast<sgpp::datadriven::FitterConfigurationDensityEstimation&>(*config));
    models.push_back(std::move(model));

    // Count the number of instances for this class
    classNumberInstances.push_back(0u);

    std::cout << "Registered new class label " << label << " with index " << idx << std::endl;
    return idx;
  } else {
    return classIdx.at(label);
  }
}

MultiGridRefinementFunctor *ModelFittingClassification::getRefinementFunctor(
    std::vector<Grid*> grids, std::vector<DataVector*> surpluses) {
  sgpp::base::AdpativityConfiguration& refinementConfig = this->config->getRefinementConfig();
  switch (refinementConfig.refinementFunctorType) {
    case RefinementFunctorType::Surplus : {
      return new MultiSurplusRefinementFunctor(grids, surpluses, refinementConfig.noPoints_,
          refinementConfig.levelPenalize, refinementConfig.threshold_);
    }
    case RefinementFunctorType::ZeroCrossing : {
      return new ZeroCrossingRefinementFunctor(grids, surpluses,
        refinementConfig.noPoints_, refinementConfig.levelPenalize,
        refinementConfig.precomputeEvaluations);
    }
    case RefinementFunctorType::DataBased : {
      if (refinementConfig.scalingCoefficients.size() != 0) {
        if (refinementConfig.scalingCoefficients.size() < models.size()) {
          std::string errorMessage = "Not enough scaling coefficients were given for the amount"
              "of classes";
          throw new application_exception(errorMessage.c_str());
        } else if (refinementConfig.scalingCoefficients.size() > models.size()) {
          std::cout << "Did not train on at least one sample for every class. Data based " <<
              "refinement not possible in this iteration..." << std::endl;
          return nullptr;
        }
      }
      return new DataBasedRefinementFunctor(grids, surpluses, &(dataset->getData()),
          &(dataset->getTargets()), refinementConfig.noPoints_, refinementConfig.levelPenalize,
          refinementConfig.scalingCoefficients);
    }
    case RefinementFunctorType::SurplusVolume : {
      std::string errorMessage = "Unsupported refinement functor type SurplusVolume "
          "for classification!";
      throw new application_exception(errorMessage.c_str());
    }
    case RefinementFunctorType::GridPointBased : {
      return new GridPointBasedRefinementFunctor(grids, surpluses, refinementConfig.levelPenalize,
          refinementConfig.precomputeEvaluations, refinementConfig.threshold_);
    }
    case RefinementFunctorType::MultipleClass : {
      return new MultipleClassRefinementFunctor(grids, surpluses, refinementConfig.noPoints_,
          0, refinementConfig.threshold_);
    }
    default:
      return nullptr;
  }
}

bool ModelFittingClassification::refine() {
  sgpp::base::AdpativityConfiguration& refinementConfig = this->config->getRefinementConfig();
  if (refinementsPerformed < refinementConfig.numRefinements_) {
    // Assemble grids and alphas
    std::vector<Grid*> grids;
    std::vector<DataVector*> surpluses;
    grids.reserve(models.size());
    surpluses.reserve(models.size());
    for (size_t idx = 0; idx < models.size(); idx++) {
      grids.push_back(&(models[idx]->getGrid()));
      surpluses.push_back(&(models[idx]->getSurpluses()));
    }

    // Create a refinement functor
    MultiGridRefinementFunctor* func = getRefinementFunctor(grids, surpluses);

    // Apply refinements for all models
    if (func) {
      // Refinement for multiple class is fundamentaly different! this needs to be fixed!
      if (refinementConfig.refinementFunctorType == RefinementFunctorType::MultipleClass) {
        // The functor handles refinements for all grids
        MultipleClassRefinementFunctor *multifunc = dynamic_cast<MultipleClassRefinementFunctor*>(
            func);
        multifunc->refine();
      } else {
        // The refinements have to be triggered manually
        for (size_t idx = 0; idx < models.size(); idx++) {
          // Precompute evaluations in case of data based / zero crossing refinement
          if (refinementConfig.precomputeEvaluations && (
              refinementConfig.refinementFunctorType == RefinementFunctorType::DataBased ||
              refinementConfig.refinementFunctorType == RefinementFunctorType::ZeroCrossing ||
              refinementConfig.refinementFunctorType == RefinementFunctorType::GridPointBased)) {
            func->preComputeEvaluations();
          }
          func->setGridIndex(idx);
          // TODO(fuchgsdk): Interaction refinement
          // In case of multiple class refinement the refinement is organized by the functor
          grids[idx]->getGenerator().refine(*func);
        }
      }

      // Apply changes to all models
      for (size_t idx = 0; idx < models.size(); idx++) {
        // TODO(fuchsgdk): Coarsening for classification? Any criteria availible?
        std::list<size_t> coarsened;
        models[idx]->refine(grids[idx]->getSize(), &coarsened);
        std::cout << "Refined model for class index " << idx << " (new size : "
            << (grids[idx]->getSize()) << ")" << std::endl;
      }
      delete func;
    }
    refinementsPerformed++;
    return true;
  }
  return false;
}

void ModelFittingClassification::update(Dataset& newDataset) {
  fit(newDataset);
}

void ModelFittingClassification::resetState() { refinementsPerformed = 0; }

}  // namespace datadriven
}  // namespace sgpp





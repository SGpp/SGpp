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
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClassification.hpp>

#include <string>
#include <vector>
#include <map>

using sgpp::base::Grid;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

using sgpp::base::application_exception;

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

// TODO(lettrich): exceptions have to be thrown if not valid.
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

bool ModelFittingClassification::refine() {
  // TODO(fuchsgdk): obviously needs to be done
  return true;
}

void ModelFittingClassification::update(Dataset& newDataset) {
  fit(newDataset);
}

void ModelFittingClassification::resetState() { refinementsPerformed = 0; }

}  // namespace datadriven
}  // namespace sgpp





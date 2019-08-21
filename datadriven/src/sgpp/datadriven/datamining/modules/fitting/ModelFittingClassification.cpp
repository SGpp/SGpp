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

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClassification.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/datadriven/configuration/DensityEstimationConfiguration.hpp>
#include <sgpp/datadriven/datamining/configuration/RefinementFunctorTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCG.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCombi.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOffParallel.hpp>
#include <sgpp/datadriven/functors/MultiSurplusRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/GridPointBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/MultipleClassRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>

#include <list>
#include <map>
#include <string>
#include <vector>

#include <fstream>
#include <iostream>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;

using sgpp::base::application_exception;

using sgpp::datadriven::RefinementFunctorType;

namespace sgpp {
namespace datadriven {

ModelFittingClassification::ModelFittingClassification(
    const FitterConfigurationClassification& config)
    : refinementsPerformed{0} {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));

#ifdef USE_SCALAPACK
  auto& parallelConfig = this->config->getParallelConfig();
  if (parallelConfig.scalapackEnabled_) {
    processGrid = std::make_shared<BlacsProcessGrid>(config.getParallelConfig().processRows_,
                                                     config.getParallelConfig().processCols_);
  }
#endif
}

double ModelFittingClassification::evaluate(const DataVector& sample) {
  if (models.size() == 0) {
    std::string errorMessage = "Prediction impossible! No models were trained!";
    throw application_exception(errorMessage.c_str());
  } else {
    auto& learnerConfig = this->config->getLearnerConfig();
    double prediction = 0.0, maxDensity = 0.0;

    // Pre compute the total number of instances
    size_t numInstances = 0;
    for (auto& p : classIdx) {
      size_t idx = p.second;
      numInstances += classNumberInstances[idx];
    }

    bool evaluatedModel = false;
    for (auto& p : classIdx) {
      double label = p.first;
      size_t idx = p.second;
      if (classNumberInstances[idx] == 0) {
        // The model for this class was not trained -> no prediction possible for this model
        continue;
      }
      double classConditionalDensity = models[idx]->evaluate(sample);
      double prior;
      if (learnerConfig.usePrior) {
        // Prior is realtive frequency of instances of this class
        prior = static_cast<double>(classNumberInstances[idx]) / static_cast<double>(numInstances);
      } else {
        // Uniform prior
        prior = 1.0;
      }
      double density = prior * classConditionalDensity;

      if (!evaluatedModel || density > maxDensity) {
        maxDensity = density;
        prediction = label;
      }
      evaluatedModel = true;
    }
    return prediction;
  }
}

void ModelFittingClassification::evaluate(DataMatrix& samples, DataVector& results) {
#ifdef USE_SCALAPACK
  auto& parallelConfig = this->config->getParallelConfig();
  if (parallelConfig.scalapackEnabled_) {
    if (!processGrid->isProcessInGrid()) {
      return;
    }
    DataVectorDistributed resultsDistributed(processGrid, results.size(),
                                             parallelConfig.rowBlockSize_);
    DataVector tmp(samples.getNcols());

    for (size_t i = 0; i < resultsDistributed.getLocalRows(); i++) {
      size_t globalRow = resultsDistributed.localToGlobalRowIndex(i);
      samples.getRow(globalRow, tmp);
      resultsDistributed.getLocalPointer()[i] = evaluate(tmp);
    }

    // only gather to master, as master calculates score and broadcasts it to workers
    resultsDistributed.toLocalDataVector(results);
    return;
  }
#endif  // USE_SCALAPACK

#pragma omp parallel for
  for (size_t i = 0; i < samples.getNrows(); i++) {
    DataVector tmp(samples.getNcols());
    samples.getRow(i, tmp);
    results.set(i, evaluate(tmp));
  }
}

void ModelFittingClassification::fit(Dataset& newDataset) {
  reset();
  update(newDataset);
}

std::unique_ptr<ModelFittingDensityEstimation> ModelFittingClassification::createNewModel(
    sgpp::datadriven::FitterConfigurationDensityEstimation& densityEstimationConfig) {
  if (densityEstimationConfig.getGridConfig().generalType_ ==
      base::GeneralGridType::ComponentGrid) {
    return std::make_unique<ModelFittingDensityEstimationCombi>(densityEstimationConfig);
  }
  switch (densityEstimationConfig.getDensityEstimationConfig().type_) {
    case DensityEstimationType::CG: {
      return std::make_unique<ModelFittingDensityEstimationCG>(densityEstimationConfig);
    }
    case DensityEstimationType::Decomposition: {
#ifdef USE_SCALAPACK
      if (densityEstimationConfig.getParallelConfig().scalapackEnabled_) {
        return std::make_unique<ModelFittingDensityEstimationOnOffParallel>(densityEstimationConfig,
                                                                            processGrid);
      }
#endif  // USE_SCALAPACK
      return std::make_unique<ModelFittingDensityEstimationOnOff>(densityEstimationConfig);
    }
  }
}

size_t ModelFittingClassification::labelToIdx(double label) {
  if (classIdx.find(label) == classIdx.end()) {
    // New class
    size_t idx = classIdx.size();
    classIdx[label] = idx;

    // Create a new model
    std::unique_ptr<ModelFittingDensityEstimation> model = createNewModel(
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

MultiGridRefinementFunctor* ModelFittingClassification::getRefinementFunctor(
    std::vector<Grid*> grids, std::vector<DataVector*> surpluses) {
  sgpp::base::AdaptivityConfiguration& refinementConfig = this->config->getRefinementConfig();
  switch (refinementConfig.refinementFunctorType) {
    case RefinementFunctorType::Surplus: {
      return new MultiSurplusRefinementFunctor(grids, surpluses, refinementConfig.noPoints_,
                                               refinementConfig.levelPenalize,
                                               refinementConfig.threshold_);
    }
    case RefinementFunctorType::ZeroCrossing: {
      return new ZeroCrossingRefinementFunctor(grids, surpluses, refinementConfig.noPoints_,
                                               refinementConfig.levelPenalize,
                                               refinementConfig.precomputeEvaluations);
    }
    case RefinementFunctorType::DataBased: {
      if (refinementConfig.scalingCoefficients.size() != 0) {
        if (refinementConfig.scalingCoefficients.size() < models.size()) {
          std::string errorMessage =
              "Not enough scaling coefficients were given for the amount"
              "of classes";
          throw new application_exception(errorMessage.c_str());
        } else if (refinementConfig.scalingCoefficients.size() > models.size()) {
          std::cout << "Did not train on at least one sample for every class. Data based "
                    << "refinement not possible in this iteration..." << std::endl;
          return nullptr;
        }
      }
      return new DataBasedRefinementFunctor(grids, surpluses, &(dataset->getData()),
                                            &(dataset->getTargets()), refinementConfig.noPoints_,
                                            refinementConfig.levelPenalize,
                                            refinementConfig.scalingCoefficients);
    }
    case RefinementFunctorType::SurplusVolume: {
      std::string errorMessage =
          "Unsupported refinement functor type SurplusVolume "
          "for classification!";
      throw new application_exception(errorMessage.c_str());
    }
    case RefinementFunctorType::GridPointBased: {
      return new GridPointBasedRefinementFunctor(grids, surpluses, refinementConfig.noPoints_,
                                                 refinementConfig.levelPenalize,
                                                 refinementConfig.precomputeEvaluations,
                                                 refinementConfig.threshold_);
    }
    case RefinementFunctorType::MultipleClass: {
      return new MultipleClassRefinementFunctor(grids, surpluses, refinementConfig.noPoints_, 0,
                                                refinementConfig.threshold_);
    }
  }
}

bool ModelFittingClassification::refine() {
  if (config->getGridConfig().generalType_ == base::GeneralGridType::ComponentGrid) {
    for (size_t i = 0; i < models.size(); i++) {
      models.at(i)->refine();
    }
    refinementsPerformed++;
    return true;
  }
  sgpp::base::AdaptivityConfiguration& refinementConfig = this->config->getRefinementConfig();
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
        MultipleClassRefinementFunctor* multifunc =
            dynamic_cast<MultipleClassRefinementFunctor*>(func);
        multifunc->refine();
      } else {
        // The refinements have to be triggered manually
        for (size_t idx = 0; idx < models.size(); idx++) {
          // Precompute evaluations in case of data based / zero crossing refinement
          if (refinementConfig.precomputeEvaluations &&
              (refinementConfig.refinementFunctorType == RefinementFunctorType::DataBased ||
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
        std::cout << "Refined model for class index " << idx
                  << " (new size : " << (grids[idx]->getSize()) << ")" << std::endl;
      }
      delete func;
    }
    refinementsPerformed++;
    return true;
  }
  return false;
}

void ModelFittingClassification::update(Dataset& newDataset) {
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
    DataMatrix* samples = p.second;
    models[idx]->update(*samples);
    classNumberInstances[idx] += samples->getNrows();
    delete samples;
  }
}

void ModelFittingClassification::reset() {
  models.clear();
  classNumberInstances.clear();
  classIdx.clear();
  refinementsPerformed = 0;
}

void ModelFittingClassification::storeClassificator() {
  std::cout << "Storing Classificator..." << std::endl;

  // store labels
  std::string labels;
  for (const auto& p : classIdx) {
    labels = labels + std::to_string(p.first) + ", " + std::to_string(p.second) + "\n";
  }
  std::ofstream labelsFile;
  // add the path of your labels.txt file here, in which the labels should be stored
  std::string pathToLabelsFile = "";
  labelsFile.open(pathToLabelsFile);
  labelsFile << labels;
  labelsFile.close();

  // store instances
  std::string instances;
  for (size_t i = 0; i < classNumberInstances.size(); i++) {
    instances = instances + std::to_string(classNumberInstances[i]) + "\n";
  }
  std::ofstream instancesFile;
  // add the path of your instances.txt file here, in which the instances should be stored
  std::string pathToInstancesFile = "";
  instancesFile.open(pathToInstancesFile);
  instancesFile << instances;
  instancesFile.close();

  // store grids and alphas
  std::string classificatorFile;
  std::string classificator;
  for (size_t i = 0; i < models.size(); i++) {
    classificator = "";
    classificatorFile = "";
    // add the path of you Grid_AlphaX.txt file here, in which the grids and alphas should be stored
    std::string pathToGridAlphaFile = "";
    classificator = classificator + pathToGridAlphaFile + "Grid_Alpha" + std::to_string(i) + ".txt";
    classificatorFile = classificatorFile + models[i]->storeFitter();
    std::ofstream file;
    file.open(classificator);
    file << classificatorFile;
    file.close();
  }
}

#ifdef USE_SCALAPACK
std::shared_ptr<BlacsProcessGrid> ModelFittingClassification::getProcessGrid() const {
  return processGrid;
}
#endif

}  // namespace datadriven
}  // namespace sgpp

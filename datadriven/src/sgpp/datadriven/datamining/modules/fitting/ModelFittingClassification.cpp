// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClassification.hpp>

#include <sgpp/base/exception/application_exception.hpp>

#include <sgpp/base/grid/RefinementConfiguration.hpp>
#include <sgpp/base/grid/RefinementFunctorTypeParser.hpp>

#include <sgpp/datadriven/configuration/DensityEstimationConfiguration.hpp>

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCG.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCombi.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOffParallel.hpp>

#include <sgpp/datadriven/functors/MultiSurplusRefinementFunctor.hpp>

#include <sgpp/datadriven/functors/classification/ClassificationRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/GridPointBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/MultipleClassRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>

#include <string>
#include <fstream>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <vector>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::application_exception;
using sgpp::base::RefinementFunctorType;

namespace sgpp {
namespace datadriven {

ModelFittingClassification::ModelFittingClassification(
    const FitterConfigurationClassification& config)
    : refinementsPerformed(0), initialGridSize(0) {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));

  this->objectStore = std::make_shared<DBMatObjectStore>();
  this->hasObjectStore = true;
#ifdef USE_SCALAPACK
  auto& parallelConfig = this->config->getParallelConfig();
  if (parallelConfig.scalapackEnabled_) {
    processGrid = std::make_shared<BlacsProcessGrid>(
        config.getParallelConfig().processRows_,
        config.getParallelConfig().processCols_);
  }
#endif
}

ModelFittingClassification::ModelFittingClassification(
    const FitterConfigurationClassification& config,
    std::shared_ptr<DBMatObjectStore> objectStore)
    : ModelFittingClassification(config) {
  this->objectStore = objectStore;
  this->hasObjectStore = true;
}

double ModelFittingClassification::evaluate(const DataVector& sample) {
  if (models.size() == 0) {
    std::string errorMessage = "Prediction impossible! No models were trained!";
    throw application_exception(errorMessage.c_str());
  } else {
    double prediction = 0.0, maxDensity = 0.0;

    // Pre compute the total number of instances
    size_t numInstances = 0;
    for (auto& p : classIdx) {
      size_t idx = p.second;
      numInstances += classNumberInstances[idx];
    }

    bool evaluatedModel = false;
    std::vector<double> priors = getClassPriors();
    for (auto& p : classIdx) {
      double label = p.first;
      size_t idx = p.second;
      if (classNumberInstances[idx] == 0) {
        // The model for this class was not trained -> no prediction possible
        // for this model
        continue;
      }
      double classConditionalDensity = models[idx]->evaluate(sample);
      double density = priors[idx] * classConditionalDensity;

      if (!evaluatedModel || density > maxDensity) {
        maxDensity = density;
        prediction = label;
      }
      evaluatedModel = true;
    }
    return prediction;
  }
}

void ModelFittingClassification::evaluate(DataMatrix& samples,
                                          DataVector& results) {
  if (models.size() == 0) {
    std::string errorMessage = "Prediction impossible! No models were trained!";
    throw application_exception(errorMessage.c_str());
  }

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

    // only gather to master, as master calculates score and broadcasts it to
    // workers
    resultsDistributed.toLocalDataVector(results);
    return;
  }
#endif  // USE_SCALAPACK

  std::vector<double> priors = getClassPriors();
  std::vector<DataVector> classResults(models.size());
  for (auto& p : classIdx) {
    size_t idx = p.second;
    DataVector results(samples.getNrows());
    models[idx]->evaluate(samples, results);
    results.mult(priors[idx]);
    classResults[idx] = results;
  }
  for (size_t j = 0; j < samples.getNrows(); j++) {
    double maxDensity = std::numeric_limits<double>::lowest();
    double prediction = 0.0;
    for (auto& p : classIdx) {
      size_t idx = p.second;
      if (maxDensity < classResults[idx][j]) {
        maxDensity = classResults[idx][j];
        prediction = p.first;
      }
    }
    results.set(j, prediction);
  }
}

std::vector<double> ModelFittingClassification::getClassPriors() const {
  auto& learnerConfig = this->config->getLearnerConfig();
  size_t numInstances = 0;
  for (auto& p : classIdx) {
    size_t idx = p.second;
    numInstances += classNumberInstances[idx];
  }

  std::vector<double> priors(models.size());
  for (auto& p : classIdx) {
    size_t idx = p.second;
    if (learnerConfig.usePrior_) {
      // Prior is realtive frequency of instances of this class
      priors[idx] = static_cast<double>(classNumberInstances[idx]) /
                    static_cast<double>(numInstances);
    } else {
      // Uniform prior
      priors[idx] = 1.0;
    }
  }

  return priors;
}

void ModelFittingClassification::fit(Dataset& newDataset) {
  reset();
  update(newDataset);
}

std::unique_ptr<ModelFittingDensityEstimation>
ModelFittingClassification::createNewModel(
    sgpp::datadriven::FitterConfigurationDensityEstimation&
        densityEstimationConfig) {
  if (densityEstimationConfig.getGridConfig().generalType_ ==
      base::GeneralGridType::ComponentGrid) {
    if (this->hasObjectStore) {
      return std::make_unique<ModelFittingDensityEstimationCombi>(
          densityEstimationConfig, this->objectStore);
    } else {
      return std::make_unique<ModelFittingDensityEstimationCombi>(
          densityEstimationConfig);
    }
  }
  switch (densityEstimationConfig.getDensityEstimationConfig().type_) {
    case DensityEstimationType::CG: {
      return std::make_unique<ModelFittingDensityEstimationCG>(
          densityEstimationConfig);
    }
    case DensityEstimationType::Decomposition: {
#ifdef USE_SCALAPACK
      if (densityEstimationConfig.getParallelConfig().scalapackEnabled_) {
        return std::make_unique<ModelFittingDensityEstimationOnOffParallel>(
            densityEstimationConfig, processGrid);
      }
#endif  // USE_SCALAPACK
      return std::make_unique<ModelFittingDensityEstimationOnOff>(
          densityEstimationConfig);
    }
  }

  throw application_exception("Unknown density estimation type");
}

size_t ModelFittingClassification::labelToIdx(double label) {
  if (classIdx.find(label) == classIdx.end()) {
    // New class
    size_t idx = classIdx.size();
    classIdx[label] = idx;

    // Create a new model
    std::unique_ptr<ModelFittingDensityEstimation> model = createNewModel(
        dynamic_cast<sgpp::datadriven::FitterConfigurationDensityEstimation&>(
            *config));
    models.push_back(std::move(model));

    // Count the number of instances for this class
    classNumberInstances.push_back(0u);

    std::cout << "Registered new class label " << label << " with index " << idx
              << std::endl;
    return idx;
  } else {
    return classIdx.at(label);
  }
}

MultiGridRefinementFunctor* ModelFittingClassification::getRefinementFunctor(
    std::vector<Grid*> grids, std::vector<DataVector*> surpluses,
    std::vector<double> priors) {
  sgpp::base::AdaptivityConfiguration& refinementConfig =
      this->config->getRefinementConfig();
  switch (refinementConfig.refinementFunctorType_) {
    case RefinementFunctorType::Surplus: {
      return new MultiSurplusRefinementFunctor(
          grids, surpluses, refinementConfig.numRefinementPoints_,
          refinementConfig.levelPenalize_,
          refinementConfig.refinementThreshold_);
    }
    case RefinementFunctorType::ZeroCrossing: {
      return new ZeroCrossingRefinementFunctor(
          grids, surpluses, priors, refinementConfig.numRefinementPoints_,
          refinementConfig.levelPenalize_,
          refinementConfig.precomputeEvaluations_);
    }
    case RefinementFunctorType::DataBased: {
      if (refinementConfig.scalingCoefficients_.size() != 0) {
        if (refinementConfig.scalingCoefficients_.size() < models.size()) {
          std::string errorMessage =
              "Not enough scaling coefficients were given for the amount"
              "of classes";
          throw application_exception(errorMessage.c_str());
        } else if (refinementConfig.scalingCoefficients_.size() >
                   models.size()) {
          std::cout << "Did not train on at least one sample for every class. "
                       "Data based "
                    << "refinement not possible in this iteration..."
                    << std::endl;
          return nullptr;
        }
      }
      return new DataBasedRefinementFunctor(
          grids, surpluses, priors, &(dataset->getData()),
          &(dataset->getTargets()), refinementConfig.numRefinementPoints_,
          refinementConfig.levelPenalize_, refinementConfig.scalingCoefficients_);
    }
    case RefinementFunctorType::SurplusVolume: {
      std::string errorMessage =
          "Unsupported refinement functor type SurplusVolume "
          "for classification!";
      throw application_exception(errorMessage.c_str());
    }
    case RefinementFunctorType::GridPointBased: {
      return new GridPointBasedRefinementFunctor(
          grids, surpluses, priors, refinementConfig.numRefinementPoints_,
          refinementConfig.levelPenalize_,
          refinementConfig.precomputeEvaluations_,
          refinementConfig.refinementThreshold_);
    }
    case RefinementFunctorType::MultipleClass: {
      return new MultipleClassRefinementFunctor(
          grids, surpluses, priors, refinementConfig.numRefinementPoints_, 0,
          refinementConfig.refinementThreshold_);
    }
    case RefinementFunctorType::Classification: {
      return new ClassificationRefinementFunctor(
          grids, surpluses, priors, refinementConfig.numRefinementPoints_,
          refinementConfig.numCoarseningPoints_, true,
          refinementConfig.thresholdType_,
          refinementConfig.refinementThreshold_,
          refinementConfig.coarseningThreshold_,
          refinementConfig.coarsenInitialPoints_, this->initialGridSize);
    }
  }

  return nullptr;
}

bool ModelFittingClassification::adapt() {
  std::vector<std::vector<size_t>> deletedPoints(models.size());
  if (config->getGridConfig().generalType_ ==
      base::GeneralGridType::ComponentGrid) {
    for (size_t i = 0; i < models.size(); i++) {
      models.at(i)->adapt();
    }
    refinementsPerformed++;
    return true;
  }
  sgpp::base::AdaptivityConfiguration& refinementConfig =
      this->config->getRefinementConfig();
  if (refinementsPerformed < refinementConfig.numRefinements_) {
    // Assemble grids and alphas
    std::vector<Grid*> grids;
    std::vector<DataVector*> surpluses;
    std::vector<double> priors;
    grids.reserve(models.size());
    surpluses.reserve(models.size());
    bool usePrior = this->config->getLearnerConfig().usePrior_;
    size_t numInstances = 0;
    for (auto& p : classIdx) {
      size_t idx = p.second;
      numInstances += classNumberInstances[idx];
    }
    for (size_t idx = 0; idx < models.size(); idx++) {
      grids.push_back(&(models[idx]->getGrid()));
      surpluses.push_back(&(models[idx]->getSurpluses()));
      if (usePrior) {
        priors.push_back(static_cast<double>(classNumberInstances[idx]) /
                         static_cast<double>(numInstances));
      } else {
        priors.push_back(1.0);
      }
    }

    // Create a refinement functor
    MultiGridRefinementFunctor* func =
        getRefinementFunctor(grids, surpluses, priors);

    // Apply refinements for all models
    if (func) {
      // Refinement for multiple class is fundamentaly different! this needs to
      // be fixed!
      if (refinementConfig.refinementFunctorType_ ==
          RefinementFunctorType::MultipleClass) {
        // The functor handles refinements for all grids
        MultipleClassRefinementFunctor* multifunc =
            dynamic_cast<MultipleClassRefinementFunctor*>(func);
        multifunc->refine();
      } else if (refinementConfig.refinementFunctorType_ ==
                 RefinementFunctorType::Classification) {
        ClassificationRefinementFunctor* classfunc =
            dynamic_cast<ClassificationRefinementFunctor*>(func);
        deletedPoints = classfunc->adaptAllGrids();
      } else {
        // The refinements have to be triggered manually
        for (size_t idx = 0; idx < models.size(); idx++) {
          // Precompute evaluations in case of data based / zero crossing
          // refinement
          if (refinementConfig.precomputeEvaluations_ &&
              (refinementConfig.refinementFunctorType_ ==
                   RefinementFunctorType::DataBased ||
               refinementConfig.refinementFunctorType_ ==
                   RefinementFunctorType::ZeroCrossing ||
               refinementConfig.refinementFunctorType_ ==
                   RefinementFunctorType::GridPointBased)) {
            func->preComputeEvaluations();
          }
          func->setGridIndex(idx);
          // TODO(fuchgsdk): Interaction refinement
          // In case of multiple class refinement the refinement is organized by
          // the functor
          GeometryConfiguration geometryConfig = config->getGeometryConfig();
          if (!geometryConfig.stencils_.empty()) {
            GridFactory gridFactory;
            grids[idx]->getGenerator().refineInter(
                *func, gridFactory.getInteractions(geometryConfig));
          } else {
            grids[idx]->getGenerator().refine(*func);
          }
        }
      }

      // Apply changes to all models
      for (size_t idx = 0; idx < models.size(); idx++) {
        models[idx]->adapt(grids[idx]->getSize(), deletedPoints.at(idx));
        std::cout << "Refined model for class index " << idx
                  << " (new size : " << (grids[idx]->getSize()) << ")"
                  << std::endl;
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

  // save initial size of the grid, if not set already
  if (this->initialGridSize == 0) {
    this->initialGridSize = models[0]->getGrid().getSize();
  }
}

void ModelFittingClassification::reset() {
  models.clear();
  classNumberInstances.clear();
  classIdx.clear();
  refinementsPerformed = 0;
}

void ModelFittingClassification::resetTraining() {
  for (auto& model : models) {
    model->resetTraining();
  }
}

void ModelFittingClassification::updateRegularization(double lambda) {
  for (auto& model : models) {
    model->updateRegularization(lambda);
  }
}

void ModelFittingClassification::storeClassificator() {
  std::cout << "Storing Classificator..." << std::endl;

  // store labels
  std::string labels;
  for (const auto& p : classIdx) {
    labels = labels + std::to_string(p.first) + ", " +
             std::to_string(p.second) + "\n";
  }
  std::ofstream labelsFile;
  // add the path of your labels.txt file here, in which the labels should be
  // stored
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
  // add the path of your instances.txt file here, in which the instances should
  // be stored
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
    // add the path of you Grid_AlphaX.txt file here, in which the grids and
    // alphas should be stored
    std::string pathToGridAlphaFile = "";
    classificator = classificator + pathToGridAlphaFile + "Grid_Alpha" +
                    std::to_string(i) + ".txt";
    classificatorFile = classificatorFile + models[i]->storeFitter();
    std::ofstream file;
    file.open(classificator);
    file << classificatorFile;
    file.close();
  }
}

std::vector<std::unique_ptr<ModelFittingDensityEstimation>>*
ModelFittingClassification::getModels() {
  return &(models);
}

std::map<double, size_t> ModelFittingClassification::getClassIdx() {
  return this->classIdx;
}
#ifdef USE_SCALAPACK
std::shared_ptr<BlacsProcessGrid> ModelFittingClassification::getProcessGrid()
    const {
  return processGrid;
}
#endif

}  // namespace datadriven
}  // namespace sgpp

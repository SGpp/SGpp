// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClustering.hpp>


using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;

namespace sgpp {
namespace datadriven {

ModelFittingClustering::ModelFittingClustering(
  const FitterConfigurationClustering& config)
  : refinementsPerformed{0} {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationClustering>(config));

#ifdef USE_SCALAPACK
    auto& parallelConfig = this->config->getParallelConfig();
  if (parallelConfig.scalapackEnabled_) {
    processGrid = std::make_shared<BlacsProcessGrid>(config.getParallelConfig().processRows_,
                                                     config.getParallelConfig().processCols_);
  }
#endif
}

void ModelFittingClustering::fit(Dataset& newDataset) {
  return;
}

void ModelFittingClustering::update(Dataset& dataset) {
  return;
}

double ModelFittingClustering::evaluate(const DataVector &sample) {
  return 0.0;
}

void ModelFittingClustering::evaluate(DataMatrix& samples, DataVector& results) {
  return;
}

bool ModelFittingClustering::refine() {
  return true;
}

void ModelFittingClustering::reset() {
  return;
}

std::unique_ptr<ModelFittingDensityEstimation>*
    ModelFittingClustering::getDensityEstimationModel() {
  return &densityEstimationModel;
}

std::unique_ptr<ModelFittingClassification>*
    ModelFittingClustering::getClassificationModel() {
  return &classificationModel;
}

}  // namespace datadriven
}  // namespace sgpp

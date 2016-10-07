/*
 * Scorer.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Metric.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ShufflingFunctor.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <memory>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Superclass for TestError, TrainError, CrossValidation, AIC, etc.
 */
class Scorer {
 public:
  Scorer(Metric* metric, ShufflingFunctor* shuffling, int64_t seed = -1);
  virtual ~Scorer();
  virtual double calculateScore(ModelFittingBase& model, Dataset& dataset,
                                double* stdDeviation = nullptr) = 0;

 protected:
  virtual void randomizeIndices(std::vector<size_t>& randomizedIndices, size_t size);
  virtual void splitSet(Dataset& fullDataset, Dataset& trainDataset, Dataset& testDataset,
                        size_t trainSize, size_t testSize,
                        const std::vector<size_t>& randomizedIndices, size_t offset = 0);
  double train(ModelFittingBase& model, Dataset& trainDataset, Dataset& testDataset);

  std::shared_ptr<Metric> metric;
  std::shared_ptr<ShufflingFunctor> shuffling;
};

} /* namespace datadriven */
} /* namespace sgpp */

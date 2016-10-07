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
  Scorer(Metric* metric, ShufflingFunctor* shuffling, int64_t seed = -1)
      : metric(std::shared_ptr<Metric>(metric)),
        shuffling(std::shared_ptr<ShufflingFunctor>(shuffling)) {
    if (seed != -1) {
      shuffling->setSeed(seed);
    }
  };

  virtual ~Scorer(){};
  virtual double calculateScore(ModelFittingBase& model, Dataset& dataset,
                                double* stdDeviation = nullptr) = 0;

 protected:
  std::shared_ptr<Metric> metric;
  std::shared_ptr<ShufflingFunctor> shuffling;
};

} /* namespace datadriven */
} /* namespace sgpp */

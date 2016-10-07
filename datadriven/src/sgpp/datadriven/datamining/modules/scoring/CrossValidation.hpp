/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * CrossValidation.hpp
 *
 *  Created on: 31.07.2016
 *      Author: Michael Lettrich
 */

#pragma once
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Metric.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ShufflingFunctor.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <memory>

namespace sgpp {
namespace datadriven {

class CrossValidation {
 public:
  virtual ~CrossValidation();
  CrossValidation(Metric* metric, ShufflingFunctor* shuffling, int64_t seed = -1,
                  size_t foldNumber = 5);

  double calculateScore(ModelFittingBase& model, Dataset& dataset, double* stdDeviation = nullptr);

 private:
  std::shared_ptr<Metric> metric;
  std::shared_ptr<ShufflingFunctor> shuffling;
  size_t foldNumber;
};

} /* namespace datadriven */
} /* namespace sgpp */

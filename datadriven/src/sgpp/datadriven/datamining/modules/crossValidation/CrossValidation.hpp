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
#include <sgpp/datadriven/datamining/modules/crossValidation/ShufflingFunctor.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Metric.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <memory>

namespace sgpp {
namespace datadriven {

class CrossValidation {
 public:
  virtual ~CrossValidation();
  CrossValidation(std::shared_ptr<Metric> metric, std::shared_ptr<ShufflingFunctor> shuffling,
                  int64_t seed = -1);

  double calculateScore(ModelFittingBase& model, Dataset& dataset, size_t foldNumber = 5,
                        std::shared_ptr<double> stdDeviation = nullptr);

  // TODO (lettrich): implement calculateScore
  // void optimizeHyperparameters ( model, dataset, hyperParamtersFixed, hyperParametersSearched,
  // fold_number, );

 private:
  std::shared_ptr<Metric> metric;
  std::shared_ptr<ShufflingFunctor> shuffling;
};

} /* namespace datadriven */
} /* namespace sgpp */

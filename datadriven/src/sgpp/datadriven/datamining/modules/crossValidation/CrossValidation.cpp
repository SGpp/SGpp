/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * CrossValidation.cpp
 *
 *  Created on: 31.07.2016
 *      Author: Michael Lettrich
 */

#include "CrossValidation.hpp"

namespace sgpp {
namespace datadriven {

CrossValidation::~CrossValidation() {
  // TODO(Michael Lettrich): Auto-generated destructor stub
}

CrossValidation::CrossValidation(const Metric& metric, const ShufflingFunctor& shuffling,
                                 int64_t seed) {  // TODO(Michael Lettrich): implement
}

double CrossValidation::calculateScore(const ModelFittingBase& model, const Dataset& dataset,
                                       size_t foldNumber, std::shared_ptr<double> stdDeviation) {
  // TODO(Michael Lettrich): implement
  return 0;
}

} /* namespace datadriven */
} /* namespace sgpp */

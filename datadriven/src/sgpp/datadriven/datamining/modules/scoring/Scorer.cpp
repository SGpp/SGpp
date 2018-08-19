/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * Scorer.cpp
 *
 * Created on: Oct 7, 2016
 *     Author: Michael Lettrich
 */

#include "Scorer.hpp"
#include <vector>

namespace sgpp {
namespace datadriven {

Scorer::Scorer(Metric* metric)
    : metric{std::unique_ptr<Metric>{metric}} { }

double Scorer::test(ModelFittingBase& model, Dataset& testDataset) {
  DataVector predictedValues{testDataset.getNumberInstances()};
  model.evaluate(testDataset.getData(), predictedValues);
  // set score
  return metric->measure(predictedValues, testDataset.getTargets());
}

} /* namespace datadriven */
} /* namespace sgpp */

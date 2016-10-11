/* Copyright (C) 2008-today The SG++ project

 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SplittingScorer.hpp
 *
 *  Created on:	07.10.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <memory>
#include <vector>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Metric.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Scorer.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

class SplittingScorer : public Scorer {
 public:
  SplittingScorer(Metric* metric, ShufflingFunctor* shuffling, int64_t seed = -1,
                  double trainPortion = 0.8);
  virtual ~SplittingScorer();

  virtual double calculateScore(ModelFittingBase& model, Dataset& dataset,
                                double* stdDeviation = nullptr);

 private:
  double trainPortion;
};

} /* namespace datadriven */
} /* namespace sgpp */

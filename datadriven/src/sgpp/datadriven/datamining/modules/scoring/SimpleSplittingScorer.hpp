/*
 * SimpleSplittingScorer.hpp
 *
 *  Created on: Feb 8, 2016
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
#include "../../configuration/DataMiningConfiguration.hpp"

 namespace sgpp {
 namespace datadriven {

 class SimpleSplittingScorer : public Scorer {
 public:
  SimpleSplittingScorer(std::shared_ptr<Metric> metric, std::shared_ptr<ModelFittingBase> fitter,
                        datadriven::DataMiningConfiguration config);
  virtual ~SimpleSplittingScorer();

  virtual double getScore(Dataset& dataset) override;

 protected:
  void splitset(Dataset& dataset, Dataset& trainingSet, Dataset& testSet, bool permute = true);

 private:
  std::unique_ptr<Dataset> trainDataset;
  std::unique_ptr<Dataset> testDataset;
  double trainPortion;
  long int seed;
};

} /* namespace datadriven */
} /* namespace sgpp */

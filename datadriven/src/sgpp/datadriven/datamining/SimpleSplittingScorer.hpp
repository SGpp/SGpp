/*
 * SimpleSplittingScorer.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#pragma once

#include <memory>
#include <vector>

#include <sgpp/datadriven/datamining/Metric.hpp>
#include <sgpp/datadriven/datamining/ModelFittingBase.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/datamining/Scorer.hpp>
#include <sgpp/datadriven/datamining/DataMiningConfiguration.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

class SimpleSplittingScorer : public Scorer {
 public:
  SimpleSplittingScorer(std::shared_ptr<Metric> metric,
                        std::shared_ptr<ModelFittingBase> fitter,
                        datadriven::DataMiningConfiguration config);
  virtual ~SimpleSplittingScorer();

  virtual double getScore(Dataset& dataset) override;

 protected:
  void splitset(Dataset& dataset, Dataset& trainingSet, Dataset& testSet,
                bool permute = true);

 private:
  Dataset trainDataset;
  Dataset testDataset;
  double trainPortion;
  long int seed;
};

} /* namespace datadriven */
} /* namespace SGPP */
//git please do not delete me

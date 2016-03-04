/*
 * SparseGridMiner.hpp
 *
 *  Created on: Feb 9, 2016
 *      Author: franzefn
 */

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/datamining/configuration/DataMiningConfiguration.hpp>
#include <sgpp/datadriven/datamining/dataSource/SampleProvider.hpp>
#include <sgpp/datadriven/datamining/scoring/Scorer.hpp>

namespace sgpp {
namespace datadriven {

class SparseGridMiner {
 public:
  SparseGridMiner(datadriven::DataMiningConfiguration pconfig);
  virtual ~SparseGridMiner();

  void run();

 private:
  datadriven::SampleProvider* dataset;
  datadriven::Scorer* scorer;

  datadriven::DataMiningConfiguration& config;
};

} /* namespace datadriven */
} /* namespace sgpp */

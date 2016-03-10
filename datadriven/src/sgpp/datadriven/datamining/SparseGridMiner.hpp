/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
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
  explicit SparseGridMiner(datadriven::DataMiningConfigJsonParser pconfig);
  virtual ~SparseGridMiner();

  void run();

 private:
  datadriven::Scorer* scorer;
  datadriven::DataMiningConfigJsonParser& config;
  datadriven::SampleProvider* dataset;
};

} /* namespace datadriven */
} /* namespace sgpp */

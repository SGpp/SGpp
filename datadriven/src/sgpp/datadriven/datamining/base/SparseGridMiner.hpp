/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SparseGridMiner.hpp
 *
 * Created on: Oct 7, 2016
 *     Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Scorer.hpp>

#include <memory>

namespace sgpp {
namespace datadriven {

class SparseGridMiner {
 public:
  SparseGridMiner(DataSource* dataSource, ModelFittingBase* fitter, Scorer* scorer);
  virtual ~SparseGridMiner();

  void learn();

 private:
  std::unique_ptr<DataSource> dataSource;
  std::unique_ptr<ModelFittingBase> fitter;
  std::unique_ptr<Scorer> scorer;
};

} /* namespace datadriven */
} /* namespace sgpp */

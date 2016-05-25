/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SparseGridMiner.hpp
 *
 *  Created on: Feb 9, 2016
 *      Author: franzefn, Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/base/DmModule.hpp>
#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/datamining/base/DmStateStorage.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>

namespace sgpp {
namespace datadriven {

class SparseGridMiner : public DmModule {
 public:
  SparseGridMiner();
  virtual ~SparseGridMiner();
  virtual void run();

 private:
  std::unique_ptr<DataSource> dataSource;
  std::unique_ptr<ModelFittingBase> fitter;
  std::shared_ptr<DmStateStorage> sharedState;
};

} /* namespace datadriven */
} /* namespace sgpp */

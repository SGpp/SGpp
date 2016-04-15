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

namespace sgpp {
namespace datadriven {

class SparseGridMiner : public DmModule {
 public:
  SparseGridMiner();
  explicit SparseGridMiner(std::shared_ptr<DmModule> firstModule,
                           std::shared_ptr<DmStateStorage> sharedState);
  virtual ~SparseGridMiner();
  void run();
};

} /* namespace datadriven */
} /* namespace sgpp */

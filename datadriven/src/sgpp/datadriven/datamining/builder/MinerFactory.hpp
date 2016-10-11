/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * MinerFactory.hpp
 *
 * Created on: Oct 10, 2016
 *     Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>

namespace sgpp {
namespace datadriven {

class MinerFactory {
 public:
  MinerFactory(){};
  virtual ~MinerFactory(){};

  virtual SparseGridMiner* buildMiner(const std::string& path) = 0;
};

} /* namespace datadriven */
} /* namespace sgpp */

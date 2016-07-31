/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ShufflingFunctor.hpp
 *
 *  Created on: 31.07.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <memory>
#include <vector>

namespace sgpp {
namespace datadriven {

class ShufflingFunctor {
 public:
  ShufflingFunctor();
  virtual void operator()(std::shared_ptr<std::vector<size_t>>) = 0;
  virtual ~ShufflingFunctor();
};

} /* namespace datadriven */
} /* namespace sgpp */

/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * Module.hpp
 *
 *  Created on: Apr 12, 2016
 *      Author: Michael Lettrich
 */

#pragma once

namespace sgpp {
namespace datadriven {

class DmModule {
 public:
  virtual ~DmModule() {}
  virtual void run() = 0;
};
} /* namespace datadriven */
} /* namespace sgpp */

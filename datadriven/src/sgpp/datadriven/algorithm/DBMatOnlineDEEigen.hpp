/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOnlineDEEigen.hpp
 *
 *  Created on: Apr 8, 2017
 *      Author: Michael Lettrich
 */

#ifdef USE_GSL

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>

namespace sgpp {
namespace datadriven {

using sgpp::datadriven::DataVector;

/**
 * TODO(lettrich) : write documentation
 */

class DBMatOnlineDEEigen : public DBMatOnlineDE {
 public:
  explicit DBMatOnlineDEEigen(DBMatOffline& offline, double beta = 0.);

 protected:
  void solveSLE(DataVector& b, bool do_cv) override;
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* USE_GSL */

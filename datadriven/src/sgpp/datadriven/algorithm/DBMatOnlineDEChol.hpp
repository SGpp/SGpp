/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOnlineDEChol.hpp
 *
 *  Created on: Apr 8, 2017
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>

#include <sgpp/datadriven/algorithm/DBMatDMSChol.hpp>

/**
 * Class that stores, generates and manipulates a density function during online phase in on/off
 * learning. This specialization operates on offline objects based on different Cholesky
 * decompositions.
 */
namespace sgpp {
namespace datadriven {

using sgpp::datadriven::DataVector;

class DBMatOnlineDEChol : public DBMatOnlineDE {
 public:
  explicit DBMatOnlineDEChol(DBMatOffline& offline, double beta = 0.);

 protected:
  void solveSLE(DataVector& b, bool do_cv) override;

  DBMatDMSChol* buildCholSolver(DBMatOffline& offlineObject, bool doCV) const;
};

} /* namespace datadriven */
} /* namespace sgpp */

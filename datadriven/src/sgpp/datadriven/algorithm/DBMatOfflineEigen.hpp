/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOfflineEigen.hpp
 *
 *  Created on: 02.03.2017
 *      Author: Michael Lettrich
 */
#ifdef USE_GSL

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>

#include <gsl/gsl_permutation.h>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * DBMatOffline specialization that uses a eigen factorization on
 * a dense matrix. Eigen factorization allows cheap changing of the regularization parameter.
 */
class DBMatOfflineEigen : public DBMatOffline {
 public:
  explicit DBMatOfflineEigen(const DBMatDensityConfiguration& oc);

  explicit DBMatOfflineEigen(const std::string& fileName);

  DBMatOffline* clone() override;

  bool isRefineable() override;

  void decomposeMatrix() override;
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* USE_GSL */

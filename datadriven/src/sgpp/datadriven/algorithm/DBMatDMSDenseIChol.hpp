/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatDMSDenseIChol.hpp
 *
 *  Created on: Apr 16, 2017
 *      Author: Michael Lettrich
 */

#pragma once
#include <sgpp/datadriven/algorithm/DBMatDMSChol.hpp>

#include <sgpp/base/grid/Grid.hpp>

namespace sgpp {
namespace datadriven {

using sgpp::base::Grid;
using sgpp::base::DataMatrix;

class DBMatDMSDenseIChol : public DBMatDMSChol {
 public:
  DBMatDMSDenseIChol(Grid* grid, double lambda, bool doCV);

 protected:
  void choleskyUpdateLambda(sgpp::base::DataMatrix& decompMatrix, double lambdaUp) const override;

 private:
  void updateProxyMatrixLambda(double lambda_up) const;

  Grid* grid;
  mutable DataMatrix proxyMatrix;
};

} /* namespace datadriven */
} /* namespace sgpp */

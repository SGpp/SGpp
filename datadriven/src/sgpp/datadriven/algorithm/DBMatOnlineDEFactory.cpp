/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOnlineDEFactory.cpp
 *
 *  Created on: Apr 8, 2017
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp>

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEEigen.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDELU.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE_SMW.hpp>

namespace sgpp {
namespace datadriven {

using sgpp::base::factory_exception;

DBMatOnlineDE* DBMatOnlineDEFactory::buildDBMatOnlineDE(DBMatOffline& offline, Grid& grid,
                                                        double lambda, double beta,
                                                        MatrixDecompositionType matDecompType) {
  // auto decompositionType = offline.getDecompositionType();
  auto decompositionType = matDecompType;
  switch (decompositionType) {
    case MatrixDecompositionType::Eigen:
#ifdef USE_GSL
      return new DBMatOnlineDEEigen(offline, grid, lambda, beta);
#else
      throw factory_exception("built without GSL");
#endif /*USE_GSL*/

    case MatrixDecompositionType::LU:
#ifdef USE_GSL
      return new DBMatOnlineDELU(offline, grid, lambda, beta);
#else
      throw factory_exception("built without GSL");
#endif /*USE_GSL*/

    case MatrixDecompositionType::Chol:
    case MatrixDecompositionType::DenseIchol:
      return new DBMatOnlineDEChol(offline, grid, lambda, beta);

    case MatrixDecompositionType::OrthoAdapt:
#ifdef USE_GSL
      return new DBMatOnlineDEOrthoAdapt(offline, grid, lambda, beta);
#else
      throw factory_exception("built without GSL");
#endif

    case MatrixDecompositionType::SMW_chol:
    case MatrixDecompositionType::SMW_ortho:
#ifdef USE_GSL
      return new DBMatOnlineDE_SMW(offline, grid, lambda, beta);
#else
      throw factory_exception("built without GSL");
#endif
  }
}

} /* namespace datadriven */
} /* namespace sgpp */

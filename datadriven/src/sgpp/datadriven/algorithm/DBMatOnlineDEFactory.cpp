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
#include <sgpp/datadriven/algorithm/DBMatOnlineDEChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEEigen.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDELU.hpp>

namespace sgpp {
namespace datadriven {

using sgpp::base::factory_exception;

DBMatOnlineDE* DBMatOnlineDEFactory::buildDBMatOnlineDE(DBMatOffline& offline, double beta) {
  auto& config = offline.getConfig();

  switch (config.decomp_type_) {
    case DBMatDecompostionType::Eigen:
#ifdef USE_GSL
      return new DBMatOnlineDEEigen(offline, beta);
#else
      throw factory_exception("built withot GSL");
#endif /*USE_GSL*/
      break;
    case DBMatDecompostionType::LU:
#ifdef USE_GSL
      return new DBMatOnlineDELU(offline, beta);
#else
      throw factory_exception("built withot GSL");
#endif /*USE_GSL*/
      break;
    case DBMatDecompostionType::Chol:
    case DBMatDecompostionType::DenseIchol:
      return new DBMatOnlineDEChol(offline, beta);
      break;
    default:
      throw factory_exception{"Unknown decomposition type."};
  }
}

} /* namespace datadriven */
} /* namespace sgpp */

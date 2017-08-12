/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatOnlineDEOrthoAdapt.hpp
 */

#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

class DBMatOnlineDEOrthoAdapt : public DBMatOnlineDE {
 public:
  /**
   * Constructor
   * Builds DBMatOnlineDEOrthoAdapt object from given offline object
   *
   * @param offline The offline object
   */
  explicit DBMatOfflineOrthoAdapt(DBMatOffline& offline, double beta = 0.);
  // yeah
}
}  // namespace datadriven
}  // namespace sgpp

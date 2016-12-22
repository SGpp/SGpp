/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * FitterTypeParser.hpp
 *
 *  Created on: 22.12.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>

#include <map>
#include <string>

namespace sgpp {
namespace datadriven {

class FitterTypeParser {
 public:
  static FitterType parse(const std::string& input);

  static const std::string& toString(FitterType type);

 private:
  typedef std::map<FitterType, std::string> FitterTypeMap_t;

  static const FitterTypeMap_t fitterTypeMap;
};

} /* namespace datadriven */
} /* namespace sgpp */

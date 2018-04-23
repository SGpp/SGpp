
/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * GeneralGridTypeParser.hpp
 *
 * Created on: 23.04.18
 *     Author: Dominik Fuchsgruber
 */

#ifndef DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_CONFIGURATION_GENERALGRIDTYPEPARSER_HPP_
#define DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_CONFIGURATION_GENERALGRIDTYPEPARSER_HPP_


#pragma once

#include <sgpp/base/grid/Grid.hpp>

#include <map>
#include <string>

namespace sgpp {
namespace datadriven {

class GeneralGridTypeParser {
 public:
  /**
   * Parses an input string and returns the corresponding sgpp::base::GeneralGridType type.
   * Throws an exception if the string has no representation.
   */
  static sgpp::base::GeneralGridType parse(const std::string& input);

  /**
   * Returns the string representation of a sgpp::base::GeneralGridType type.
   */
  static const std::string& toString(sgpp::base::GeneralGridType type);

 private:
  typedef std::map<sgpp::base::GeneralGridType, std::string> GeneralGridTypeMap_t;

  /**
   * Map to describe the mapping between sgpp::base::GeneralGridType types and their string
   * representations.
   */
  static const GeneralGridTypeMap_t generalGridTypeMap;
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_CONFIGURATION_GENERALGRIDTYPEPARSER_HPP_ */

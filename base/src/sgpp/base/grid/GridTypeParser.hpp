// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>

#include <map>
#include <string>

namespace sgpp {
namespace base {

class GridTypeParser {
 public:
  /**
   * Convert strings to values #sgpp::base::GridType. Throws if there is no
   * valid
   * representation
   * @param input case insensitive string representation of a
   * #sgpp::base::GridType.
   * @return the corresponding #sgpp::base::GridType.
   */
  static GridType parse(const std::string &input);

  /**
   * generate string representations for values of #sgpp::base::GridType.
   * @param type enum value.
   * @return string representation of a #sgpp::base::GridType.
   */
  static const std::string &toString(GridType type);

 private:
  typedef std::map<GridType, std::string> GridTypeMap_t;

  /**
   * Map containing all values of  #sgpp::base::GridType and the corresponding
   * string representation.
   */
  static const GridTypeMap_t gridTypeMap;
};
} /* namespace base */
} /* namespace sgpp */

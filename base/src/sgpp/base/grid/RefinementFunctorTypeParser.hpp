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

class RefinementFunctorTypeParser {
 public:
  /**
   * Convert strings to values #sgpp::base::RefinementFunctorType. Throws if
   * there is no valid
   * representation
   * @param input case insensitive string representation of a
   * #sgpp::base::RefinementFunctorType.
   * @return the corresponding #sgpp::base::RefinementFunctorType.
   */
  static RefinementFunctorType parse(const std::string& input);

  /**
   * generate string representations for values of
   * #sgpp::base::RefinementFunctorType.
   * @param type enum value.
   * @return string representation of a #sgpp::base::RefinementFunctorType.
   */
  static const std::string& toString(RefinementFunctorType type);

 private:
  typedef std::map<RefinementFunctorType, std::string>
      RefinementFunctorTypeMap_t;

  /**
   * Map containing all values of #sgpp::base::RefinementFunctorType and the
   * corresponding
   * string representation.
   */
  static const RefinementFunctorTypeMap_t refinementFunctorTypeMap;
};

} /* namespace base */
} /* namespace sgpp */

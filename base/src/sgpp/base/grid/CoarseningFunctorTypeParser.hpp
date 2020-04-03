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

class CoarseningFunctorTypeParser {
 public:
  /**
   * Convert strings to values #sgpp::base::CoarseningFunctorType. Throws if there is no valid
   * representation
   * @param input case insensitive string representation of a
   * #sgpp::base::CoarseningFunctorType.
   * @return the corresponding #sgpp::base::CoarseningFunctorType.
   */
  static CoarseningFunctorType parse(const std::string& input);

  /**
   * generate string representations for values of
   * #sgpp::base::CoarseningFunctorType.
   * @param type enum value.
   * @return string representation of a #sgpp::base::CoarseningFunctorType.
   */
  static const std::string& toString(CoarseningFunctorType type);

 private:
  typedef std::map<CoarseningFunctorType, std::string> CoarseningFunctorTypeMap_t;

  /**
   * Map containing all values of #sgpp::base::CoarseningFunctorType and the corresponding string
   * representation.
   */
  static const CoarseningFunctorTypeMap_t coarseningFunctorTypeMap;
};

} /* namespace base */
} /* namespace sgpp */

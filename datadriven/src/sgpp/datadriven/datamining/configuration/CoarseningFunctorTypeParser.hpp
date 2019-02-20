#pragma once

#include <sgpp/base/grid/Grid.hpp>

#include <map>
#include <string>

namespace sgpp {
namespace datadriven {
using sgpp::base::CoarseningFunctorType;

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
   * generate string representations for values of #sgpp::base::CoarseningFunctorType.
   * @param type enum value.
   * @return string representation of a #sgpp::base::CoarseningFunctorType.
   */
  static const std::string& toString(CoarseningFunctorType type);

 private:
  typedef std::map<CoarseningFunctorType, std::string> CoarseningFunctorTypeMap_t;

  /**
   * Map containing all values of #sgpp::base::CoarseningFunctorType and the corresponding
   * string representation.
   */
  static const CoarseningFunctorTypeMap_t coarseningFunctorTypeMap;
};

} /* namespace datadriven */
} /* namespace sgpp */

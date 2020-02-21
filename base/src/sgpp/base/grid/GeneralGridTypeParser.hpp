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

/**
 * Parser class to parse a general grid type into a GeneralGridType enum type
 * and vice versa
 */
class GeneralGridTypeParser {
 public:
  /**
   * Parses an input string and returns the corresponding
   * sgpp::base::GeneralGridType type.
   * Throws an exception if the string has no representation
   * @param input the grid type to parse
   * @return the parsed grid type
   */
  static sgpp::base::GeneralGridType parse(const std::string& input);

  /**
   * Returns the string representation of a sgpp::base::GeneralGridType type
   * @param type the grid type to retrieve the string representation from
   * @return the string representation of the the grid type
   */
  static const std::string& toString(sgpp::base::GeneralGridType type);

 private:
  // Define a type for the map grid type -> string
  typedef std::map<sgpp::base::GeneralGridType, std::string>
      GeneralGridTypeMap_t;

  /**
   * Map to describe the mapping between sgpp::base::GeneralGridType types and
   * their string
   * representations.
   */
  static const GeneralGridTypeMap_t generalGridTypeMap;
};

} /* namespace base */
} /* namespace sgpp */

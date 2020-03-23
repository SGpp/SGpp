// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/solver/TypesSolver.hpp>

#include <map>
#include <string>

namespace sgpp {
namespace solver {

/**
 * Convenience class to convert strings to #sgpp::solver::SLESolverType and generate
 * string representations for values of #sgpp::solver::SLESolverType.
 */
class SLESolverTypeParser {
 public:
  /**
   * Convert strings to values #sgpp::solver::SLESolverType. Throws if there is no valid
   * representation
   * @param input case insensitive string representation of a
   * #sgpp::solver::SLESolverType.
   * @return the corresponding #sgpp::solver::SLESolverType.
   */
  static SLESolverType parse(const std::string &input);

  /**
   * generate string representations for values of #sgpp::solver::SLESolverType.
   * @param type enum value.
   * @return string representation of a #sgpp::solver::SLESolverType.
   */
  static const std::string &toString(SLESolverType type);

 private:
  typedef std::map<SLESolverType, std::string> SLESolverTypeMap_t;

  /**
   * Map containing all values of  #sgpp::solver::SLESolverType and the corresponding
   * string representation.
   */
  static const SLESolverTypeMap_t sleSolverTypeMap;
};
} /* namespace solver */
} /* namespace sgpp */

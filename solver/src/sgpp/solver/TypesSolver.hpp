// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef TYPESSOLVER_HPP
#define TYPESSOLVER_HPP

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cstddef>

namespace sgpp {

namespace solver {

/**
 * enum to address different SLE solvers in a standardized way
 */
enum class SLESolverType { CG, BiCGSTAB };

struct SLESolverConfiguration {
  sgpp::solver::SLESolverType type_;
  double eps_;
  size_t maxIterations_;
  double threshold_;
};

struct SLESolverSPConfiguration {
  sgpp::solver::SLESolverType type_;
  float eps_;
  size_t maxIterations_;
  float threshold_;
};

class SolverTypeParser {
 public:
  static SLESolverType parse(const std::string& input) {
    auto inputLower = input;
    std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

    if (inputLower.compare("cg") == 0) {
      return sgpp::solver::SLESolverType::CG;
    } else if (inputLower.compare("bicgstab") == 0) {
      return sgpp::solver::SLESolverType::BiCGSTAB;
    } else {
      std::string errorMsg =
          "Failed to convert string \"" + input + "\" to any known SLESolverType";
      throw base::data_exception(errorMsg.c_str());
    }
  };
};
}  // namespace solver
}  // namespace sgpp

#endif /* TYPESSOLVER_HPP */

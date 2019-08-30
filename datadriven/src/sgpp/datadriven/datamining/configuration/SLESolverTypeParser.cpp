/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SolverTypeParser.cpp
 *
 * Created on: Jan 25, 2017
 *     Author: Michael Lettrich
 */

#include <sgpp/datadriven/datamining/configuration/SLESolverTypeParser.hpp>

#include <sgpp/base/exception/data_exception.hpp>
#include <algorithm>
#include <string>

namespace sgpp {
namespace datadriven {

using sgpp::solver::SLESolverType;

SLESolverType sgpp::datadriven::SLESolverTypeParser::parse(const std::string &input) {
  auto inputLower = input;
  std::transform(inputLower.begin(), inputLower.end(), inputLower.begin(), ::tolower);

  if (inputLower.compare("cg") == 0) {
    return sgpp::solver::SLESolverType::CG;
  } else if (inputLower.compare("bicgstab") == 0) {
    return sgpp::solver::SLESolverType::BiCGSTAB;
  } else if (inputLower.compare("fista") == 0) {
    return sgpp::solver::SLESolverType::FISTA;
  } else {
    std::string errorMsg = "Failed to convert string \"" + input + "\" to any known SLESolverType";
    throw base::data_exception(errorMsg.c_str());
  }
}

const std::string &sgpp::datadriven::SLESolverTypeParser::toString(SLESolverType type) {
  return sleSolverTypeMap.at(type);
}

const SLESolverTypeParser::SLESolverTypeMap_t SLESolverTypeParser::sleSolverTypeMap = []() {
  return SLESolverTypeParser::SLESolverTypeMap_t{std::make_pair(SLESolverType::CG, "CG"),
                                                 std::make_pair(SLESolverType::BiCGSTAB,
                                                                "BiCGSTAB"),
                                                 std::make_pair(SLESolverType::FISTA, "FISTA")};
}();
} /* namespace datadriven */
} /* namespace sgpp */
